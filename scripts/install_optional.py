# scripts/install_optional.py

import subprocess

# Optional packages to install via pip
optional_packages = [
    "git+https://github.com/bp-kelley/descriptastorus",
    "meeko",
    "prody"
]

for package in optional_packages:
    print(f"🔧 Installing: {package}")
    subprocess.run(["pip", "install", package], check=True)

print("✅ Optional packages installed successfully.")
