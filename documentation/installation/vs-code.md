# Visual Studio Code

You will need to install a source code editor so you are able to work with the Process code. Visual 
Studio Code, or VSCode for short, is a great choice. It is a lightweight but powerful source
code editor which runs on your desktop and is available for Windows, MacOS and Linux. It has a
vast extension package allowing ease of use with a range of languages. Visual Studio Code is
available for install [here](https://code.visualstudio.com/) and more information can be found on
the [docs pages](https://code.visualstudio.com/docs).

-------------

## WSL with VS Code

To connect WSL to VS Code, the following extension should be installed in VS Code upon
opening: `Remote - WSL`. This allows for opening of any folder in the Windows Subsystem for
Linux. This is performed by using `Ctrl+Shift+X` on VS Code, searching for `Remote - WSL` and
then installing.

--------------

## Automatically activating virtual environment on VS Code Open

When VS Code is first opened, you are able to set it such that the command:

```bash
source env/bin/activate
```

is executed automatically. This saves manually activating the virtual environment every time you 
open the application. This is done by first using `Ctrl+Shift+P` and searching for 
`Python:Select Interpreter`. The select: `Python *version* ('env':venv) ./.venv/bin/python`. This 
should be starred as the recommended version.

You will see that in your project the .vscode directory will contain a `settings.json` file. Open 
this and inside of it add:

```json
"python.terminal.activateEnvironment": true
```

Don't forget to add a comma before to separate it from already present key value pairs.

Now, close your terminal and close VS Code. Reopen and open a new terminal which should now 
automatically point to the virtual environment signalled by an `(.venv)` in front of your user.
