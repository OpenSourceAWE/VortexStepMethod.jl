## Notes for developers

### How to check the .cff file

You can validate the `CITATION.cff` file with the command:
```bash
cffconvert --validate
```

You can install the script `cffconvert` with the command:
```bash
python3 -m venv venv
source venv/bin/activate
python3 -m pip install cffconvert
```
After using the tool, deactivate the Python environment with:
```bash
deactivate
```