"""Entry point shared by the installed command and source-directory launcher."""


def main(argv=None):
    from .cli import main as cli_main
    return cli_main(argv)
