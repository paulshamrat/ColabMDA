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

## Code Quality Standards

We use the following tools to maintain code quality:
*   **Black:** For code formatting.
*   **Ruff:** For linting and static analysis.

Happy coding!
