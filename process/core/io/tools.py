import importlib

import click


def mfile_opt(exists: bool = False):
    return click.option(
        "-f",
        "--mfile",
        "mfile",
        default="MFILE.DAT",
        type=click.Path(exists=exists),
        help="The mfile to read",
    )


mfile_arg = click.argument("mfiles", nargs=-1, type=click.Path(exists=True))


def indat_opt(default="IN.DAT"):
    return click.option(
        "-i",
        "--input",
        "indat",
        type=click.Path(exists=True),
        help="The path to the input file",
        default=default,
    )


def save(help_):
    return click.option("-s", "--save", "save", default=False, is_flag=True, help=help_)


def split_callback(ctx: click.Context, param, value: str | None) -> list[str] | None:  # noqa: ARG001
    return value.split(":") if isinstance(value, str) else value


### Taken from click documentation
class LazyGroup(click.Group):
    def __init__(self, *args, lazy_subcommands=None, **kwargs):
        super().__init__(*args, **kwargs)
        # lazy_subcommands is a map of the form:
        #
        #   {command-name} -> {module-name}.{command-object-name}
        #
        self.lazy_subcommands = lazy_subcommands or {}

    def list_commands(self, ctx):
        base = super().list_commands(ctx)
        lazy = sorted(self.lazy_subcommands.keys())
        return sorted(base + lazy)

    def get_command(self, ctx, cmd_name):
        if cmd_name in self.lazy_subcommands:
            return self._lazy_load(cmd_name)
        return super().get_command(ctx, cmd_name)

    def _lazy_load(self, cmd_name):
        # lazily loading a command, first get the module name and attribute name
        import_path = self.lazy_subcommands[cmd_name]
        modname, cmd_object_name = import_path.rsplit(".", 1)
        # do the import
        mod = importlib.import_module(modname)
        # get the Command object from that module
        cmd_object = getattr(mod, cmd_object_name)
        # check the result to make debugging easier
        if not isinstance(cmd_object, click.Command):
            raise ValueError(
                f"Lazy loading of {import_path} failed by returning a non-command object"
            )
        return cmd_object


###
