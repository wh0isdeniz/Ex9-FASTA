import random
import csv

def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000) -> int:
    """Gets an integer from the user within [min_val, max_val].
    Repeats the prompt on invalid input instead of crashing."""
    while True:
        raw = input(prompt)
        try:
            value = int(raw)
            if min_val <= value <= max_val:
                return value
            else:
                print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")
        except ValueError:
            print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")


