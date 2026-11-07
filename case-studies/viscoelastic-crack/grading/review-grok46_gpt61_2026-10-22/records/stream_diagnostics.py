"""Transparent Unix relay retaining SSE types/error categories, never bodies."""
from __future__ import annotations

import argparse
from collections import Counter
import datetime as dt
import json
import os
from pathlib import Path
import select
import socket
import socketserver
import threading
import time


class Events:
    def __init__(self):
        self.buffer = b''
        self.types = Counter()
        self.errors = []

    def feed(self, data):
        self.buffer += data
        while b'\n' in self.buffer:
            line, self.buffer = self.buffer.split(b'\n', 1)
            if not line.startswith(b'data:'):
                continue
            try:
                event = json.loads(line[5:].strip())
            except (ValueError, UnicodeError):
                continue
            if not isinstance(event, dict):
                continue
            kind = str(event.get('type', 'unknown'))[:100]
            self.types[kind] += 1
            response = event.get('response', {})
            errors = [event.get('error')]
            if isinstance(response, dict):
                errors.append(response.get('error'))
            for error in errors:
                if not isinstance(error, dict):
                    continue
                message = str(error.get('message', '')).lower()
                safe = {k: str(error[k])[:100] for k in ('code', 'type', 'param')
                        if isinstance(error.get(k), (str, int))}
                safe['message_terms'] = [term for term in (
                    'timeout', 'context', 'token', 'limit', 'reasoning', 'encrypted',
                    'signature', 'history', 'tool', 'function', 'duplicate',
                    'required', 'invalid', 'unsupported', 'internal', 'server',
                    'rate', 'filter', 'content', 'empty', 'stream', 'id',
                ) if term in message]
                safe['event_type'] = kind
                if safe not in self.errors:
                    self.errors.append(safe)


class Relay(socketserver.BaseRequestHandler):
    def handle(self):
        started = time.monotonic()
        events = Events()
        upstream = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        upstream.connect(str(self.server.upstream))
        record = {'started_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
                  'request_bytes': 0, 'response_bytes': 0}
        headers = b''
        in_body = False
        try:
            while True:
                ready, _, _ = select.select([self.request, upstream], [], [], 3600)
                if not ready:
                    record['stop_reason'] = 'idle limit'
                    return
                for source in ready:
                    data = source.recv(65536)
                    if not data:
                        record['stop_reason'] = 'upstream EOF' if source is upstream else 'client EOF'
                        return
                    if source is upstream:
                        record['response_bytes'] += len(data)
                        if in_body:
                            events.feed(data)
                        else:
                            headers += data
                            if b'\r\n\r\n' in headers:
                                head, body = headers.split(b'\r\n\r\n', 1)
                                record['http_status_line'] = head.split(b'\r\n', 1)[0].decode('ascii', 'replace')[:100]
                                headers = b''
                                in_body = True
                                events.feed(body)
                    else:
                        record['request_bytes'] += len(data)
                    (self.request if source is upstream else upstream).sendall(data)
        except (ConnectionError, OSError) as error:
            record['transport_error_type'] = type(error).__name__
        finally:
            upstream.close()
            record.update(elapsed_seconds=time.monotonic()-started,
                          event_types=dict(events.types), errors=events.errors,
                          content_retained=False)
            with self.server.log.open('a') as out:
                out.write(json.dumps(record)+'\n')


class Server(socketserver.ThreadingUnixStreamServer):
    daemon_threads = True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('socket_dir', type=Path)
    parser.add_argument('run_json', type=Path)
    parser.add_argument('log', type=Path)
    args = parser.parse_args()
    active = args.socket_dir/'model.sock'
    original = args.socket_dir/'model-original.sock'
    temporary = args.socket_dir/'model-diagnostic.sock'
    if original.exists() or temporary.exists():
        raise RuntimeError('A diagnostic relay already exists')
    server = Server(str(temporary), Relay)
    server.upstream, server.log = original, args.log
    os.chmod(temporary, 0o666)
    active.rename(original)
    temporary.rename(active)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        while json.loads(args.run_json.read_text()).get('status') == 'running':
            time.sleep(5)
    finally:
        server.shutdown()
        server.server_close()


if __name__ == '__main__':
    main()
