from logging import Handler

import process.process_output as process_output
from process.fortran import constants


class ProcessLogHandler(Handler):
    def __init__(self, capturing=True):
        super().__init__()
        self._logs = []
        self._capturing = capturing

    def start_capturing(self):
        self._capturing = True

    def stop_capturing(self):
        self._capturing = False

    def emit(self, record):
        if self._capturing:
            self._logs.append(record)

    def clear_logs(self):
        self._logs = []

    def num_logs(self):
        return len(self._logs)

    def render_warnings(self):
        return "\n\n".join([
            f"({w.pathname}:{w.lineno}) {w.getMessage()}" for w in self._logs
        ])


logging_model_handler = ProcessLogHandler()


def show_errors(file_unit: int):
    warning_string = (
        "******************************************** Errors and Warnings *********************************************"
        f"\n{logging_model_handler.render_warnings()}"
    )
    print(warning_string)
    process_output.write(file_unit, warning_string)
    process_output.ovarre(
        constants.mfile,
        "Error status",
        "(error_status)",
        0 if logging_model_handler.num_logs() == 0 else 2,
        "OP",
    )
