import os, sys, pty, timeit, subprocess

# These methods are all really close in timing.
# Looks like it doesn't matter which we use for speed.

COMMAND = ['ls', '-l', '/Users/alapoint/Documents/CAMRA']

class PseudoTerminalPair():
    
    def __init__(self, command: list[str]):
        self.controller, self.worker = pty.openpty()
        self.run_process(command)
    
    def run_process(self, command):
        try:
            process = subprocess.Popen(command, stdin=self.worker, 
                            stdout=self.worker, stderr=self.worker, 
                            shell=False)
            process.wait()
        except Exception as e:
            print(f"There was an error running the command '{command}':\n{e}")
    
    def get_controller(self):
        return self.controller
    
    def get_worker(self):
        return self.worker
    
    def get_output(self):
        output = os.read(self.controller, 1024).decode()
        return output
    
    def close_terminal(self):
        os.close(self.controller)
        os.close(self.worker)


def with_pty_openpty():
    controller, worker = pty.openpty()
    process = subprocess.Popen(COMMAND, stdin=worker, stdout=worker, stderr=worker, shell=False)
    process.wait()
    output = os.read(controller, 1024).decode()
    os.close(controller)
    os.close(worker)


def with_pty_object():
    process = PseudoTerminalPair(COMMAND)
    output = process.get_output()
    process.close_terminal()


def with_subprocess_run():
    process = subprocess.run(COMMAND, capture_output=True, shell=False)


def measure_runtime(func, iterations=100):
    timer = timeit.Timer(func)
    total_time = timer.timeit(number=iterations)
    return total_time / iterations


def main():
    iterations = 1000

    avg_time1 = measure_runtime(with_pty_openpty, iterations)
    avg_time2 = measure_runtime(with_pty_object, iterations)
    avg_time3 = measure_runtime(with_subprocess_run, iterations)

    print(f"Average runtime of with_pty_openpty over {iterations} iterations: {avg_time1:.6f} seconds")
    print(f"Average runtime of with_pty_object over {iterations} iterations: {avg_time2:.6f} seconds")
    print(f"Average runtime of with_subprocess_run over {iterations} iterations: {avg_time3:.6f} seconds")

if __name__ == "__main__":
    main()
