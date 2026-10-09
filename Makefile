.PHONY: setup pipeline dashboard

setup:
	@command -v uv >/dev/null 2>&1 || python3 -m pip install uv
	uv sync

pipeline:
	uv run python load_data.py
	uv run python run_analysis.py

dashboard:
	uv run streamlit run dashboard.py --server.port 8501 --server.address 0.0.0.0 --server.headless true