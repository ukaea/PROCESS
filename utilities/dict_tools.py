import importlib.util
import pathlib

CURRENT_DIR = pathlib.Path(__file__).resolve().parent


def module_from_file(module_name, file_path):
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


proc_dict = module_from_file(
    "create_dicts_config", CURRENT_DIR.parent / "scripts/create_dicts_config.py"
)
