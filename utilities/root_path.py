from pathlib import Path
import time

timeout = time.time() + 10  # 2 minutes from now
found_root = False
back = "../"
while not found_root:
    if time.time() > timeout:
        print(
            "Can't find repository root. Make sure utility is being run "
            "inside a PROCESS repository clone"
        )
        break
    else:
        my_file = Path(back + ".gitignore")
        if my_file.is_file():
            found_root = True
            REPO_ROOT = back
        back += "../"
