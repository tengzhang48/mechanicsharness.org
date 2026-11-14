"""A fresh opencode home for one agent run: a config that sets that model's thinking to its maximum and denies the web
tools, the package cache, and beside it `HOME_DIR.env` with the key of ONE provider as environment variables for
`docker run --env-file` (opencode 2.2 keeps logins in its database, not in auth.json, so the key is handed over
through the environment). Nothing else of the organizer's home. Same as the magnet case's script, plus the web rule.

    python3 make_opencode_home.py HOME_DIR azure|deepseek|zai-coding-plan
Models (13 Nov 2026): azure/Kimi-K3.1-Pro (reasoningEffort max); deepseek/deepseek-flash = "DeepSeek V4.1 Flash"
(variant max: thinking enabled, reasoningEffort max); zai-coding-plan/glm-5.3 (opencode's own variant max).
Web rule (Teng, 13 Nov 2026: "Web tools off, log audited"): the actions webfetch and websearch are denied for every
agent of opencode, sub-agents included; --auto does not approve a denied action. The network itself stays on, because
the model is reached through it, so each run's log is audited afterwards for any other use of the network.
The .env file is never printed and must be deleted after the launch.
"""
import json, os, shutil, sys
home, provider = sys.argv[1], sys.argv[2]
src = os.path.expanduser('~')
auth = json.load(open(os.path.join(src, '.local/share/opencode/auth.json')))
os.makedirs(os.path.join(home, '.config/opencode'))
a = auth[provider]
env = {'azure': lambda: {'AZURE_RESOURCE_NAME': a['metadata']['resourceName'], 'AZURE_API_KEY': a['key']},
       'deepseek': lambda: {'DEEPSEEK_API_KEY': a['key']}, 'zai-coding-plan': lambda: {'ZHIPU_API_KEY': a['key']}}[provider]()
p = home.rstrip('/') + '.env'
fd = os.open(p, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
with os.fdopen(fd, 'w') as f:
    f.write(''.join('%s=%s\n' % kv for kv in env.items()))
MODELS = {
    "azure": {"models": {"Kimi-K3.1-Pro": {"name": "Kimi K3.1 Pro (Azure Kimi-K3.1-Pro, reasoning max)",
                                           "options": {"reasoningEffort": "max", "useCompletionUrls": True},
                                           "variants": {"max": {"reasoningEffort": "max"}}}}},
    "deepseek": {"models": {"deepseek-flash": {"reasoning": True, "options": {"thinking": {"type": "enabled"}, "reasoningEffort": "max"},
                                               "variants": {"max": {"thinking": {"type": "enabled"}, "reasoningEffort": "max"}}}}},
    "zai-coding-plan": None}
config = {"$schema": "https://opencode.ai/config.json", "autoupdate": False, "share": "disabled",
          "permission": {"webfetch": "deny", "websearch": "deny"}}
if MODELS[provider]:
    config["provider"] = {provider: MODELS[provider]}
json.dump(config, open(os.path.join(home, '.config/opencode/opencode.jsonc'), 'w'), indent=1)
if os.path.isdir(os.path.join(src, '.cache/opencode')):
    shutil.copytree(os.path.join(src, '.cache/opencode'), os.path.join(home, '.cache/opencode'), symlinks=True)
print('opencode home for', provider, 'in', home)
