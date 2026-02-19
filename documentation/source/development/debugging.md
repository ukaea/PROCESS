# Debugging with PDB

Debugging allows you stop a program mid-execution (called 'breaking') and explore the current state (e.g. print the current value of variables or attributes). In Python, this is done using `pdb` which comes as-standard with Python. 

First, you must decide where in the code you would like to break. When debugging, a useful way to 
decide is to look at the terminal for where the error may have occurred and break there. To insert a 
breakpoint in the code, insert the following code where you want the code to break:

```python
import pdb
pdb.set_trace()
```

This will add a breakpoint so that the code will stop running on the line after this import
statement. You are able to insert multiple break statements in the code if you choose. Next, you 
should rerun your previous command for which you experienced an error (although debugging can be used whenever, not just when an error occurs). The code will run to the 
breakpoint and then break in the terminal. Effectively, the running of the code has 'paused' such 
that you are able to look around in your current state. This means that in the terminal you are able 
to print out the values of variables at the current moment to understand if they are correct or to 
diagnose if one is giving an error. This is especially useful when looking for 0 or infinity 
errors. From here you can decide what to do next.

A print statement would be along the following lines:

```bash
(pdb): print(variablename)
```

### What else can you do when you are in the break?

So you have broken at the break point, what can you do now other than print? There are two key 
commands that are useful: 'c' and 's' 

- 'c' is continue, and when run at a breakpoint this will continue the execution beyond the current 
  breakpoint to the next or to the end of the program.
- 's' is step. This executes the current line, stops at the first possible occasion- whether that 
  is in a function that is called or the next line in the current function.

Other interesting commands can be found in Python's PDB documentation 
[here](https://docs.python.org/3/library/pdb.html). Essentially, you are able to write any Python 
code in the interactive terminal.

To exit the debugger just use:

```bash
(pdb): quit() or q
```

### Using the VS Code GUI to Debug

You are also able to use VS Code's built in visual debugger to debug your code rather than using 
break statements inserted in the code itself. For an in depth instruction guide on how to do this, 
see instructions [here](https://code.visualstudio.com/docs/editor/debugging).
