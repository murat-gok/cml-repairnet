.PHONY: setup doctor data stage1 labels features train stage25-quick stage25-full stage3-quick stage3-full stage4-doctor stage4-quick stage4-full stage4-audit test lint all virtual-help

setup:
	bash scripts/bootstrap.sh

doctor:
	cml-repairnet doctor

data:
	cml-repairnet download

stage1:
	cml-repairnet stage1 --no-download

labels:
	python scripts/prepare_bsa_labels.py

features:
	python scripts/build_sequence_features.py
	python scripts/build_structure_features.py

train:
	python scripts/train_baselines.py

stage25-quick:
	python -u scripts/run_stage25.py --mode quick

stage25-full:
	python -u scripts/run_stage25.py --mode full --bootstrap 20000 --permutations 50000 --permutation-repeats 100 --n-jobs 1

stage3-quick:
	python scripts/run_stage3.py --mode quick --n-jobs 1

stage3-full:
	python scripts/run_stage3.py --mode full --n-jobs 1

stage4-doctor:
	python scripts/run_stage4.py doctor

stage4-quick:
	python -u scripts/run_stage4.py all --mode quick --workers 2

stage4-full:
	python -u scripts/run_stage4.py all --mode full --workers 4 --top-n 100

stage4-audit:
	python scripts/run_stage4.py audit

virtual-help:
	python scripts/select_virtual_simulation_batch.py --help
	python scripts/run_external_protein_validation.py --help

test:
	python -m pytest -q

lint:
	ruff check src scripts tests

all: data stage1 train
