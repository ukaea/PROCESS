"""The module containing PROCESS' custom log handler.

This handler is designed to capture logs during the output phase of
a run--when the models are being run for the final time.
This ensures captured logs are relevant for the solution that is written
to the MFile.
"""

from logging import Handler

import process.process_output as process_output
from process import constants


class ProcessLogHandler(Handler):
    """PROCESS' custom log handler that stores captured logs on the handler.

    The handler will only store logs when an internal attribute _capturing is
    True. This can be changed by setting the capturing keyword when instantiating
    the handler or using the methods start_capturing/stop_capturing.
    """

    def __init__(self, capturing=True):
        """Instantiates a ProcessLogHandler.

        :param capturing: capture and store emitted logs?
        :type capturning: bool
        """
        super().__init__()
        self._logs = []
        self._capturing = capturing

    def start_capturing(self):
        """Start capturing logs and storing them within the handler."""
        self._capturing = True

    def stop_capturing(self):
        """Stop capturing logs and storing them within the handler."""
        self._capturing = False

    def emit(self, record):
        """Method called when creating a logging record (e.g. logger.warning)."""
        if self._capturing:
            self._logs.append(record)

    def clear_logs(self):
        """Empties the handler's internal record of previous logs."""
        self._logs = []

    def num_logs(self):
        """The number of logs the handler has stored."""
        return len(self._logs)

    def render_warnings(self):
        """Render the stored warnings for printing to the terminal or OUTFile."""
        return "\n\n".join([
            f"({w.pathname}:{w.lineno}) {w.getMessage()}" for w in self._logs
        ])


logging_model_handler = ProcessLogHandler()


def show_errors(file_unit: int):
    """Write the rendered captured logs to the terminal/OUTFile

    Parameters
    ----------
    file_unit : int
        a number describing the output medium (terminal, OUTFile)

    """
    warning_string = (
        "******************************************** Errors and Warnings *********************************************"
        f"\n{logging_model_handler.render_warnings()}"
    )
    print(warning_string)
    process_output.write(file_unit, warning_string)
    process_output.ovarre(
        constants.MFILE,
        "Error status",
        "(error_status)",
        0 if logging_model_handler.num_logs() == 0 else 2,
        "OP",
    )
