# Frequently Asked Questions

## I have installed XLaunch to run a GUI and now applications including PROCESS is running very slowly

*This is intended for users who have installed WSL on a Windows machine*

Some programs require the use of XLaunch. If you have installed this, you should have followed these [instructions](https://intranet.ukaea.uk/software/guides/wsl2.html) with the section `Enable x-forwarding`. There is proved to be the issue that upon starting your laptop, if you then do NOT activate XLaunch before running PROCESS, it will be incredibly slow. This can be tested by running:

```python
pytest -k baseline
```

This should take less than or around 5 seconds. Users have experienced this taking around 135 seconds following XLaunch install.

**There are two solutions to this issue**

1. If you plan on using the XLaunch features often, then it will be easiest to simply activate XLaunch before running PROCESS. You can do this by simply running the program and following the steps.
1. If you no longer use or only infrequently use Xlaunch features, it is recommended you do the following to disable Xlaunch:
    - Go to your .bashrc file in an Ubuntu terminal using
    ```bash
    nano ~/.bashrc
    ```
    - You will see the following line near the bottom, which was inserted when downloading and installing Xlaunch:
    ```bash
    export DISPLAY=$(awk '/nameserver / {print $2; exit}' /etc/resolv.conf 2>/dev/null):0
    ```
    - Remove this line and exit the .bashrc file making sure you save the changes.
    - Restart the terminal and test that the changes have been implemented by again running the below in the PROCESS directory:
    ```python
    pytest -k baseline
    ```
    - As mentioned, this should only take a few seconds
