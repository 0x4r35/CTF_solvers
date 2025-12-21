from pwn import *
import re

context.log_level = 'info'

def solve():
    HOST = '35.194.98.181'
    PORT = 10961

    try:
        r = remote(HOST, PORT)

        r.recvuntil(b"input message (hex): ")

        payload_hex = b"Get Flag.".hex()
        r.sendline(payload_hex.encode())

        r.recvuntil(b"input signature (int): ")
        r.sendline(b"0")

        response = r.recvall(timeout=2).decode()

        if "TSGCTF" in response:
            flag = re.search(r"(TSGCTF\{.*?\})", response).group(1)
            print(flag)
        else:
            print(response)

    except Exception as e:
        print(e)
    finally:
        r.close()

if __name__ == "__main__":
    solve()

  #FLAG : TSGCTF{+he|23_was_4_f4rm3|2_had_a_|)o6_a|^_||)_8i|^_|6o_wa5_his_nam3}