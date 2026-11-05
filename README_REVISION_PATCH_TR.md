# CML-RepairNet Hakem Revizyon Test Yaması v0.8.0

Kurulum ve ilk çalıştırma:

```bash
cd /mnt/c/Users/user/RepairNet
unzip -o /mnt/c/Users/user/Downloads/CML_RepairNet_Reviewer_Revision_Tests_v0.8.0.zip -d .
source .venv/bin/activate

bash scripts/backup_before_revision.sh
python -m pip install -e ".[dev,experiment]"
python -m pytest -q
python scripts/run_revision_tests.py doctor
python -u scripts/run_revision_tests.py all \
  --mode quick \
  --models extra_trees ridge lightgbm \
  --groups 4 \
  --permutations 200 \
  --n-jobs 4 \
  2>&1 | tee results/revision_tests/revision_quick_console.txt
```

`PASS: pre-revision backup created and verified` görülmeden analiz komutlarına geçmeyin.
Ayrıntılı açıklama: `docs/REVISION_TESTS_GUIDE_TR.md`.
