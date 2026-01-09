def solve():
    # 1. Target string from Ruby source
    # "UF\\JC" means a single backslash byte (ASCII 92)
    target_str = r"Coufhlj@bixm|UF\JCjP^P<"
    target_bytes = [ord(c) for c in target_str]

    # 2. Generate the "Generator23" sequence
    # The sequence includes 2, 3, and then every number coprime to 2 and 3.
    key_sequence = [2, 3]
    num = 4
    while len(key_sequence) < 23:
        # Generator23 logic: allow numbers NOT divisible by 2 or 3
        if num % 2 != 0 and num % 3 != 0:
            key_sequence.append(num)
        num += 1

    # 3. Decrypt
    flag = []
    for t, k in zip(target_bytes, key_sequence):
        flag.append(chr(t ^ k))

    print(f"Key Sequence Used: {key_sequence}")
    print("Flag:", "".join(flag))

if __name__ == "__main__":
    solve()