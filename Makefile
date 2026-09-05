.PHONY: setup pipeline dashboard test

setup:
	python -m pip install -r requirements.txt

pipeline:
	python load_data.py
	python analysis.py

dashboard:
	python -m streamlit run dashboard.py --server.address=0.0.0.0 --server.headless=true --browser.gatherUsageStats=false

test:
	python -m pytest
