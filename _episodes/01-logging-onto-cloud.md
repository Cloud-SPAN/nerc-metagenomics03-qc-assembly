---
title: "Logging onto the Cloud"
teaching: 5
exercises: 40
questions:
- How do I connect to an AWS instance?
objectives:
- Log onto to a running instance
- Log off from a running instance
keypoints:
- You can use one set of log-in credentials for many instances
- Logging off an instance is not the same as turning off an instance
---

<script language="javascript" type="text/javascript">
function set_page_view_defaults() {
    document.getElementById('div_aws_win').style.display = 'block';
    document.getElementById('div_aws_unix').style.display = 'none';
};

function change_content_by_platform(form_control){
    if (!form_control || document.getElementById(form_control).value == 'aws_win') {
        set_page_view_defaults();
    } else if (document.getElementById(form_control).value == 'aws_unix') {
        document.getElementById('div_aws_win').style.display = 'none';
        document.getElementById('div_aws_unix').style.display = 'block';
        document.getElementById('div_hpc').style.display = 'none';
        document.getElementById('div_cyverse').style.display = 'none';
    } else {
        alert("Error: Missing platform value for 'change_content_by_platform()' script!");
    }
}

window.onload = set_page_view_defaults;
</script>

## Important Note

This lesson covers how to log into, and out of, an *already running* Amazon instance.

## Background to AWS

An Amazon Web Services (AWS) instance is a **remote computer** that runs on AWS infrastructure and that is accessible from any laptop or desktop provided you have the log-in credentials.

To log in if you are attending a tutor-led Genomics workshop you will need:
- the name of your instance and a login key file
- the shell/terminal application --- **Windows users** should have already installed the *Git Bash* shell; otherwise follow the directions in the [Setup](../setup)
- the *secure shell* (`ssh`) application, which is readily available in MacOS, Linux and Windows. **Windows users** will use `ssh` through Git Bash.

As the name implies, `ssh` provides you with a secure (encrypted) way to use a remote *shell*, as simple as this (you do not have to type this yet):

 ~~~
 $ ssh -i login-key-instanceNNN.pem  csuser@instanceNNN.cloud-span.aws.york.ac.uk
 ~~~
 {: .bash}

Where `NNN` is a three-digit number giving your instance number.


## Instructions

### Create a folder for the course
To keep things tidy and easily accessible, create a folder (or directory) to keep everything related to this course: your login key file, your notes, data, etc. If you have completed the Prenomics or Genomics course, you will have already made a `cloudspan` folder. If that is the case, you can navigate to your existing folder and move onto **Downloading your login key**.

In theory you can make your Cloud-SPAN directory anywhere in your file system but we recommend making it inside your Desktop folder, to make it easy to access.

1. **Create the folder** `cloudspan` in your *Desktop*.

   Minimise all windows until you can see your desktop background. Right click and select *New*, then *Folder*. Name the folder `cloudspan`.

   You should see a folder icon appear on your desktop with the label `cloudspan`.

   Additionally, if you enter your file explorer application you should be able to click on the *Desktop* directory at the side and see the `cloudspan` folder.

2. **Write down the absolute path** to your `cloudspan` folder.

   Find out what the absolute path is using your file manager application. Right click on the folder, or in any blank space inside the folder, and select *Properties*.

   The field called *Location* will tell you the absolute path for your folder. Once you have this written down, do not lose it! Now you can find your way back to the `cloudspan` folder whenever you need to, no matter where you are in your file structure.

> ## Reminder of Absolute vs relative paths
>There are two ways of writing a file path - absolute paths and relative paths.
>
>An absolute path contains the complete list of directories needed to locate a file on your computer. This allows you to reach the file no matter where you are. The example just given (home/docs/data/doc3.txt) is an absolute path.
>
>A relative path describes the location of a file relative to your current working directory. For example, if you were already in the folder called docs, the relative path for doc3.txt would be data/doc3.txt. There is no need to give instructions to navigate a route you have already taken.
>
>If, however, you were in the folder called docs and you wanted to open one of the .exe files, you would probably use the absolute path for that file (home/programs/.exe) to get there. This is because you have not navigated any of the route yet, so you need the full ‘address’.
>

### Download your login key file

Next we will download your unique login key file from the email you received from the Cloud-SPAN team. This type of file is called a `.pem` file. It contains a certificate which allows you to communicate with the Cloud securely. Without the `.pem` file you cannot access the Cloud.

For now we will use the file explorer to move the `.pem` file around.

1. **Find out where downloads are saved** on your computer.

   How you so this will depend on which browser you use. You can find instructions for changing your default download location in [Chrome](https://support.google.com/chrome/answer/95759?hl=en-GB&co=GENIE.Platform%3DDesktop), [Edge](https://support.microsoft.com/en-us/microsoft-edge/find-where-your-browser-is-saving-downloads-d3e83af6-68bb-aa90-3167-eeb657013902) or [Safari](https://support.apple.com/en-gb/guide/safari/sfri40598/mac).

   If you already know which folder your downloads go to, then you can skip this step.

2. **Download your login key file** to the folder you just created.

   Click on the link embedded in the email you received from the Cloud-SPAN team.

   **Mac users** may need to Click on 'download' when the file says it can't be opened.

   If your browser asks you "where do you want to download the file?", choose the `cloudspan` directory.

   Otherwise, once downloading is finished, copy and paste/drag and drop your login key file from wherever it was downloaded to your `cloudspan` folder.

### Open a Terminal and change the access permissions of your login key file

1. **Open the *cloudspan* folder you created for the course**

    Open your file manager and navigate to the `cloudspan` folder (hint: we recommended you make the folder in your *Desktop* directory - but you might have made it somewhere else). If you cannot find the folder, you can remind yourself where it is stored by looking at the absolute path you wrote down in the previous episode.

    The folder should contain the login key file we downloaded in the previous episode and nothing else.

2. **Right-click and open your machine's command line interface**

    Now we can open the command line.

    For Windows users:
    - Right click anywhere inside the blank space of the file manager, then select **Git Bash Here**.

    For Mac users:

    You have two options.

    EITHER
    - Open **Terminal** in one window and type `cd` followed by a space. Do not press enter! Now open **Finder** in another window. Drag and drop the `cloudspan` folder from the Finder to the Terminal. You should see the file path leading to your `cloudspan` folder appear. Now press enter to navigate to the folder.

    OR
    - Open **Terminal** and type `cd` followed by the absolute path that leads to your `cloudspan` folder. Press enter.  

    A new window will open - this is your command line interface, also known as the shell or the terminal. Once the terminal opens, it will display/output the **command prompt** to signal that it is ready to accept commands (instructions). The **command prompt** is 1 or 2 lines depending on your operating system (Windows, Linux, MacOS) and will be similar to the following.

    Typical command prompt for Windows Git Bash users:

    ~~~
    username@machineid MINGW64 ~
    $
    ~~~
    {: .output}

    Obviously "username" and "machineid" in the Output box above will be different when you open a terminal and will correspond to the actual username and the name of the machine you are using.

    The character `$` is the typical ending of user prompts (the ending of admin users prompts is typically `#`). Commands you type will follow the `$`.


    Typical command prompt for Linux users:

    ~~~
    username@machineid:~ $
    ~~~
    {: .output}

    Typical command prompt for MacOS users:

    ~~~
    machineid:~ username $
    ~~~
    {: .output}

    In the rest of this course we will show only the `$` to represent the prompt.

3. **Check that you are in the right folder**

    The terminal should have automatically set our `cloudspan` folder as the current working directory. This is because we asked the terminal to open from a specific location.

    You can check if the working directory is set correctly by looking at the file path which is defined to the left of your command prompt (`$`). It should display the second half of the absolute path we wrote down previously, usually starting after your computer's username, and always ending in `/cloudspan`.

    You can also check by typing the letters `ls` after the command prompt and pressing enter. This will list all the files in the working directory AKA all files in the `cloudspan` folder. In our case, this should be just one file, the login key ending in `.pem`.

4. **Change the access permissions of your login key file**

    Enter the following command to change the access permissions of your file but **replace** NN with the actual number in your file name:
    ~~~
    $ chmod 400 login-key-instanceNNN.pem 
    ~~~
    {: .bash}

    The command `chmod` (change access mode) makes your login key file accessible to you only (and non-accessible to any other potential users of your computer), a condition that is required and checked by the program `ssh` that you will use next to login to your AWS instance. You will learn about file access permissions later in the course.


If you were skipping the steps above having already made your `cloudspan` folder, here is where you shoud pay attention again.

## Login into your instance with ssh

1. Copy and paste the command in the Code box below to your *terminal*, but **replace** `NNN` with the number in your login key file name.

    Windows Git Bash users only:
    - **copy** the command the usual Windows way: (1) highlight it with the mouse pointer while pressing the mouse left button and (2) press Ctrl-v (keys Ctrl and v simultaneously).
    - but **paste** it the Linux/Unix way: by pressing the mouse middle button while hovering the mouse pointer over the Git Bash window. (Try with the right button if the middle button doesn't work.)

    ~~~
    $ ssh -i login-key-instanceNNN.pem  csuser@instanceNNN.cloud-span.aws.york.ac.uk
    ~~~
    {: .bash}

    *Be sure to replace* `NNN` *twice.* You can use the left and right arrow keys to move to where NN is.

    The `-i` option tells `ssh` the identity file containing the key to send to your AWS instance for it to check that you have access permissions to connect as an *ssh client*.

    **NB**: `ssh` looks for the identity file in the directory where you are, Desktop/cloudspan (that you moved to with the command `cd` above), which is the same directory wherein you downloaded your login key file.  


2. The terminal will display a security message, after you enter the `ssh` command, *similar* to the one below:

    ~~~
    The authenticity of host 'instanceNNN.cloud-span.aws.york.ac.uk (52.211.132.120)' can't be established.ECDSA key fingerprint is SHA256:8N054prkkCeM4GCDSsa0AUnSQw5ngBQHbOR40FqfqLg.
    Are you sure you want to continue connecting (yes/no/[fingerprint])? 
    ~~~
    {: .output}

    Type **yes** to continue and get connected to your AWS instance.

    The terminal will display a few messages and at the end the **prompt** of the **shell running** on your Linux AWS instance:

    ~~~
    ...
    csuser@instanceNNN-gc:~ $
    ~~~
    {: .output}

    Note that you did not need to give a password to login to your instance --- you are using your login-key file for authentication.

## Logging off your cloud instance

Logging off your instance is a lot like logging out of your local computer but it doesn't shut the computer off. **Be aware that AWS instances accrue charges whenever they are running, even if you are logged off**. Today, however, you do not need to worry about this!

To log off, use the `exit` command in the same terminal you connected with. This will close the connection, and your terminal will go back to showing your local computer prompt, for example:

~~~
csuser@instance05-gc:~ $ exit
~~~
{: .bash}
~~~
logout
Connection to instanceNNN.cloud-span.aws.york.ac.uk closed.
$
~~~
{: .output}

## Subsequent logins to your AWS instance

To login back to your instance, open a terminal, move to the directory you created for the course, and `ssh` as before:

~~~
$ cd Desktop/cloudspan
$ ssh -i login-key-instanceNNN.pem  csuser@instanceNNN.cloud-span.aws.york.ac.uk
~~~
{: .bash}
