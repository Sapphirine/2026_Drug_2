"""
Colab initialization helper.

Usage at the top of every notebook:

    from src.colab_init import setup
    PROJECT_ROOT = setup(repo_url="https://github.com/Sapphirine/2026_Drug_2.git")

This does:
  1. Mounts Google Drive at /content/drive
  2. Clones (or pulls) the repo into /content/drive/MyDrive/egfr_diffusion
  3. Creates all result subdirectories
  4. Installs missing dependencies
  5. Returns the absolute project root path
"""

from __future__ import annotations
import os
import sys
import subprocess


def _run(cmd: list[str], check: bool = True):
    print("$", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout[-500:])
    if check and result.returncode != 0:
        print("STDERR:", result.stderr[-500:])
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return result


def setup(
    repo_url: str,
    project_dir: str = "egfr_diffusion",
    install_deps: bool = True,
) -> str:
    """Set up Colab environment. Returns project root path."""

    # 1. Mount Google Drive (skip if not on Colab)
    try:
        from google.colab import drive
        drive.mount("/content/drive", force_remount=False)
        drive_root = "/content/drive/MyDrive"
    except ImportError:
        print("Not running on Colab — using local paths.")
        return os.path.abspath(os.path.join(os.getcwd(), ".."))

    # 2. Clone or pull repo
    project_root = os.path.join(drive_root, project_dir)
    if not os.path.exists(os.path.join(project_root, ".git")):
        os.makedirs(drive_root, exist_ok=True)
        _run(["git", "clone", repo_url, project_root])
    else:
        os.chdir(project_root)
        _run(["git", "pull"], check=False)

    os.chdir(project_root)
    if project_root not in sys.path:
        sys.path.insert(0, project_root)

    # 3. Ensure all directories exist
    for d in [
        "data/pdb", "data/pockets", "data/baselines",
        "results/generated", "results/vina_scores", "results/figures",
        "external",
    ]:
        os.makedirs(os.path.join(project_root, d), exist_ok=True)

    # 4. Install dependencies
    if install_deps:
        req_file = os.path.join(project_root, "requirements.txt")
        if os.path.exists(req_file):
            print("Installing requirements (first run ~3-5 min)...")
            subprocess.run(
                ["pip", "install", "-q", "-r", req_file],
                check=False,
            )
        # Open Babel (for Vina PDBQT conversion)
        subprocess.run(["apt-get", "install", "-y", "-qq", "openbabel"], check=False)

    print()
    print(f"PROJECT_ROOT = {project_root}")
    return project_root


def commit_and_push(message: str, paths: list[str] | None = None):
    """Push small result files (CSVs, figures) back to GitHub from Colab."""
    if paths is None:
        paths = ["results/vina_scores/*.csv", "results/figures/*.png",
                 "results/summary_table.csv"]
    _run(["git", "add"] + paths, check=False)
    _run(["git", "commit", "-m", message], check=False)
    _run(["git", "push"], check=False)
