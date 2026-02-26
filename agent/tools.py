import subprocess
import os


def run_tests():
    print("Running tests...")
    result = subprocess.run(
        ["make", "test"],
        capture_output=True,
        text=True
    )
    return result


def read_source_files(root_dirs=("src", "include")):
    collected_code = ""

    for root in root_dirs:
        if not os.path.exists(root):
            continue

        for subdir, _, files in os.walk(root):
            for file in files:
                if file.endswith(".cpp") or file.endswith(".h"):
                    path = os.path.join(subdir, file)
                    with open(path, "r", encoding="utf-8") as f:
                        collected_code += f"\n// FILE: {path}\n"
                        collected_code += f.read()
                        collected_code += "\n"

    return collected_code


def apply_patch(patch_text: str):
    patch_file = "agent_patch.diff"

    with open(patch_file, "w", encoding="utf-8") as f:
        f.write(patch_text)

    check = subprocess.run(
        ["git", "apply", "--check", patch_file],
        capture_output=True,
        text=True
    )

    if check.returncode != 0:
        print("Patch validation failed.")
        print(check.stderr)
        return False

    apply = subprocess.run(
        ["git", "apply", patch_file],
        capture_output=True,
        text=True
    )

    if apply.returncode != 0:
        print("Patch application failed.")
        print(apply.stderr)
        return False

    print("Patch applied successfully.")
    return True
