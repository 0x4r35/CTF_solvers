from pwn import *

# Set up the connection
host = '34.170.146.252'
port = 6218

def solve():
    # Connect to the remote server
    io = remote(host, port)

    try:
        # The server sends: "here is my 🍅: {p}"
        # We read until the colon and space to extract p
        io.recvuntil(b': ')
        p_output = io.recvline().strip()
        p = int(p_output)
        
        log.info(f"Received Prime p: {p}")

        # Apply Fermat's Little Theorem: a^(p-1) = 1 (mod p)
        # We set our choice to p - 1
        payload = str(p - 1)

        # The server asks: "what is your 🍅> "
        io.sendlineafter(b'> ', payload.encode())
        log.success(f"Sent payload (p-1)")

        # Retrieve the flag
        # The server checks the pow(a, choice, p) and prints the flag
        io.recvuntil(b"here is the flag: ")
        flag = io.recvline().decode().strip()
        
        print("\n" + "="*50)
        print(f"🚩 FLAG: {flag}")
        print("="*50 + "\n")
        
    except Exception as e:
        log.error(f"An error occurred: {e}")
    finally:
        io.close()

if __name__ == "__main__":
    solve()