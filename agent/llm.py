import os
import requests
from dotenv import load_dotenv

load_dotenv()

API_KEY = os.getenv("OPENROUTER_API_KEY")

if not API_KEY:
    raise ValueError("OPENROUTER_API_KEY not found in .env file")


def call_llm(prompt: str) -> str:
    url = "https://openrouter.ai/api/v1/chat/completions"

    headers = {
        "Authorization": f"Bearer {API_KEY}",
        "Content-Type": "application/json"
    }

    payload = {
        "model": "deepseek/deepseek-chat",
        "messages": [
            {
                "role": "system",
                "content": (
                    "You are a senior C++ engineer. "
                    "You generate minimal unified diff patches only. "
                    "Never include explanations. "
                    "Output raw git diff format only."
                )
            },
            {
                "role": "user",
                "content": prompt
            }
        ],
        "temperature": 0.1
    }

    response = requests.post(url, headers=headers, json=payload)

    if response.status_code != 200:
        raise Exception(f"LLM Error: {response.text}")

    return response.json()["choices"][0]["message"]["content"]
