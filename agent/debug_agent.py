from llm import call_llm
from tools import run_tests, read_source_files, apply_patch

MAX_ITERATIONS = 5


def build_prompt(test_output: str, source_code: str) -> str:
    return f"""
The following C++ tests are failing:

================ TEST OUTPUT ================
{test_output}
=============================================

Here is the project source code:

================ SOURCE CODE ================
{source_code}
=============================================

Task:
- Identify the root cause of the failure.
- Generate a minimal unified diff patch.
- Modify only necessary lines.
- Do NOT rewrite entire files.
- Output ONLY raw git diff format.
"""


def main():
    print("Starting AI Debug Agent...")

    for iteration in range(1, MAX_ITERATIONS + 1):
        print(f"\n=== Iteration {iteration}/{MAX_ITERATIONS} ===")

        result = run_tests()

        if result.returncode == 0:
            print("\nAll tests passed.")
            return

        print("\nTests failed. Sending to LLM...")

        source_code = read_source_files()

        prompt = build_prompt(
            result.stdout + "\n" + result.stderr,
            source_code
        )

        try:
            patch = call_llm(prompt)
        except Exception as e:
            print("LLM call failed:", e)
            return

        if not patch.strip():
            print("LLM returned empty patch.")
            return

        print("\nApplying patch...")
        success = apply_patch(patch)

        if not success:
            print("Patch rejected. Stopping.")
            return

    print("\nMaximum iterations reached. Manual review required.")


if __name__ == "__main__":
    main()
