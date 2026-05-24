# Contributing to ColabMDA

Thank you for your interest in contributing! We welcome help in making Molecular Dynamics more accessible and robust for everyone.

## How to Contribute

1.  **Report Bugs:** If you find a bug, please open an issue with a clear description and steps to reproduce.
2.  **Suggest Features:** Have an idea? Open an issue to discuss it!
3.  **Submit Pull Requests:**
    *   Fork the repository.
    *   Create a new branch for your feature or fix.
    *   Ensure your code follows the **Black** code style.
    *   Ensure all tests and linting pass (`ruff check .`).
    *   Submit a PR with a clear description of your changes.

## Development Setup

To set up a development environment locally:

```bash
git clone https://github.com/paulshamrat/ColabMDA.git
cd ColabMDA
mamba env create -f environment.yml
conda activate colabmda
```

## Project Architecture & Modularity Guide

ColabMDA is designed with a **Modular Command-Line Architecture**. This separates the "User Interface" (CLI) from the "Scientific Logic" (Engines).

### 1. The Core Structure
*   **`src/colabmda/cli.py`**: The central entry point. It handles argument parsing and basic path resolution. It **should not** contain heavy scientific logic.
*   **`src/colabmda/modeller/`**: Contains everything related to protein modeling.
    *   `engine.py`: The core MODELLER logic.
    *   `commands.py`: Bridge between the CLI and the engine.
*   **`src/colabmda/openmm/`**: Contains everything related to OpenMM simulations.
    *   `engine/`: Contains core simulation scripts.
    *   `modular/`: Contains individual EM, NVT, NPT, and MD steps.
    *   `commands.py`: Bridge between the CLI and the engine.

### 2. Adding a New Feature
1.  **Logic First:** Add your core logic in a new file within the appropriate engine directory (e.g., `src/colabmda/openmm/analysis_tools.py`).
2.  **Add Command Handler:** Update the `commands.py` file in that directory to expose your logic as a function.
3.  **Update CLI:** Add a new subcommand to `cli.py` that calls your command handler.
4.  **Test Locally:** Run your new command with a `test_` prefix to ensure it doesn't get tracked by Git.

### 3. Coding Style
*   **Functional Programming:** Prefer clean, standalone functions over complex classes where possible.
*   **Logging:** Use the project's logging patterns (writing to a `pipeline.log` in the output directory) to help users debug.
*   **Linting:** Always run `black .` and `ruff check . --fix` before committing.

Happy coding!
