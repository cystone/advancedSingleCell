**代码版本控制是作为程序员，必须要考虑的问题**，针对`Rstudio`可以利用`Github`进行版本控，下面对整个操作过程进行说明：

## 设定目录

在`windows`系统下，选择`Tools` --> `Global Options`，然后选择`Git/SVN`，选择`Git executable`，所以安装前提是你要有先安装`Git`，如下图所示

![image-20210329122736942](https://gitee.com/cystone2020/document/raw/master/image-20210329122736942.png)

选择git目录


 然后`Create RSA Key`

![image-20210329122752182](https://gitee.com/cystone2020/document/raw/master/image-20210329122752182.png)

![image-20210329122816918](https://gitee.com/cystone2020/document/raw/master/image-20210329122816918.png)

![image-20210329122840750](https://gitee.com/cystone2020/document/raw/master/image-20210329122840750.png)

打开`Tools`，选择`shell`，输入命令：
 `git config --global user.email "youremail@gmail.com`
 `git config --global user.name "yourname"`
 `ssh -T git@github.com`
 使用`GitHub`上的名字

![image-20210329122854117](https://gitee.com/cystone2020/document/raw/master/image-20210329122854117.png)GitHub连接

## 新建一个一个工程

新建一个工程，选择`New Directory`

![image-20210329122906337](https://gitee.com/cystone2020/document/raw/master/image-20210329122906337.png)

新建工程

然后勾选`Create a git repository`

![image-20210329122936787](https://gitee.com/cystone2020/document/raw/master/image-20210329122936787.png)

创建Git

这个时候`Rstudio`会出现`git`栏，提交到本地，只需要在`git`栏下面点击`commit`，即可提交至本地

![image-20210329122958635](https://gitee.com/cystone2020/document/raw/master/image-20210329122958635.png)GitHub提交



可以将代码保存至`GitHub`上，并且创建分支，在`GitHub`上创建一个`New respository`，命名为`test`

![image-20210329123014956](https://gitee.com/cystone2020/document/raw/master/image-20210329123014956.png)

GitHUb上创建

打开`Rstudio`中的`Shell`窗口，输入`git`命令



```csharp
git remote add origin  https://github.com/chengfeifan/test.git
git config remote.origin.url git@github.com:chengfeifan/test.git
git pull  origin master
git push  origin master
```

将`origin`重新定向
 `git remote set-url origin https://github.com/chengfeifan/test.git`

## 在本地新建一个`GitHub`上已经存在的项目

首先在`Rstudio`上新建一个`project`，选择`version control`,然后选`Clone Git Respository`，将`GitHub`上`repository`的`url`加入到选项中

![image-20210329123033548](https://gitee.com/cystone2020/document/raw/master/image-20210329123033548.png)

Paste_Image.png

然后在`shell`窗口输入
 `git config remote.origin.url git@github.com:ewenharrison/test.git`

## git中设置上游

在`git`的时候，我们会建立许多有特性的分支，建立分支的时候，如何使得远端也出现分支，需要用到下面的命令：

```bash
 git push --set-upstream origin master
```

## 常见错误和解决办法

$git push
fatal: unable to access 'https://github.com/cystone/advancedSingleCell.git/': LibreSSL SSL_connect: SSL_ERROR_SYSCALL in connection to github.com:443

```bash
git config --global --unset http.proxy
git config --global http.sslVerify false
```



