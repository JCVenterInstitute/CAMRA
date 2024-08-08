import os, sys, pty, subprocess, time


def login(user: str, password: str) -> None:
    """Logs into the BV-BRC service using the provided username and password.

    Args:
        user (str): BV-BRC Username. Defaults to sys.argv[1].
        password (str): BV-BRC Password. Defaults to sys.argv[2].
    """
    controller, worker = pty.openpty()
    login_command = ["p3-login", user]

    try:
        process = subprocess.Popen(
            login_command, stdin=worker, stdout=worker, stderr=worker, shell=False
        )

        while True:
            output = os.read(controller, 1024).decode()
            if "Password:" in output:
                # If prompted, enter password
                # Open a write handle to the controller-end of the pty to write to
                pin = os.fdopen(controller, "w")
                pin.write(
                    f"{password}\n"
                )  # NOTE: the \n is equivalent to the user hitting the "enter" key
                pin.flush()

                # Wait for login to complete
                process.wait()

                break
            # Otherwise wait for prompt
            time.sleep(0.1)

    except Exception as e:
        print(f"An error has occurred wh: {e}")

    finally:
        os.close(controller)
        os.close(worker)


login(sys.argv[1], sys.argv[2])
