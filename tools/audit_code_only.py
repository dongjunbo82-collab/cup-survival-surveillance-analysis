"""Audit tracked publication files. Prints counts, not sensitive matches."""
import pathlib,re,subprocess,sys
root=pathlib.Path(__file__).resolve().parents[1]
tracked=subprocess.check_output(['git','ls-files','-z'],cwd=root).decode().split('\0')
tracked=[x for x in tracked if x]
if not tracked:raise SystemExit('Stage the code allowlist before auditing.')
problems=[]
rules=[r'(?<![A-Za-z])[A-Za-z]:[\\/]',r'/Users/',r'/home/[^/\s]+/',
       r'(?i)\bgh[pousr]_[A-Za-z0-9_]{15,}',r'github_pat_[A-Za-z0-9_]+',
       r'-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----',
       r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}\b',
       r'(?i)\b(?:password|client_secret|api_key|access_token)\s*[:=]\s*[\"\x27][^\"\x27]{8,}']
for name in tracked:
    p=root/name
    allowed=(name in {'.gitignore','README.md','INPUT_SCHEMA.md','run_analysis.R'} or
             name.startswith('R/')and p.suffix=='.R' or
             name=='config/histology_definitions.json' or
             name.startswith(('python/','tools/'))and p.suffix=='.py' or
             name.startswith('tests/')and p.suffix=='.R')
    if not allowed:problems.append((name,'not code or permitted documentation'))
    text=p.read_text(encoding='utf-8-sig')
    # The audit's own patterns contain credential examples but no credentials.
    if name!='tools/audit_code_only.py':
        for pattern in rules:
            if re.search(pattern,text):problems.append((name,'privacy pattern matched'))
    if '\x00'in text:problems.append((name,'binary content'))
if problems:
    print('FAIL',len(problems),'issues in',len({x[0]for x in problems}),'files')
    for name,reason in problems:print(name,reason)
    sys.exit(1)
print('PASS:',len(tracked),'tracked code/configuration/documentation files; no data or document files.')
