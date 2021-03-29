**代码版本控制是作为程序员，必须要考虑的问题**，针对`Rstudio`可以利用`Github`进行版本控，下面对整个操作过程进行说明：

## 设定目录

在`windows`系统下，选择`Tools` --> `Global Options`，然后选择`Git/SVN`，选择`Git executable`，所以安装前提是你要有先安装`Git`，如下图所示

![img](https:////upload-images.jianshu.io/upload_images/1716465-ece69d10fb81a5c9.png?imageMogr2/auto-orient/strip|imageView2/2/w/456/format/webp)

选择git目录


 然后`Create RSA Key`



![img](https:////upload-images.jianshu.io/upload_images/1716465-8d63af6b7503a58b.png?imageMogr2/auto-orient/strip|imageView2/2/w/456/format/webp)

形成秘钥

![img](https:////upload-images.jianshu.io/upload_images/1716465-1519d51fc1df9b02.png?imageMogr2/auto-orient/strip|imageView2/2/w/456/format/webp)

复制公钥

![img](https:////upload-images.jianshu.io/upload_images/1716465-60edf2cf8c934e5f.png?imageMogr2/auto-orient/strip|imageView2/2/w/456/format/webp)

登陆github，添加秘钥

打开`Tools`，选择`shell`，输入命令：
 `git config --global user.email "youremail@gmail.com`
 `git config --global user.name "yourname"`
 `ssh -T git@github.com`
 使用`GitHub`上的名字

![img](https:////upload-images.jianshu.io/upload_images/1716465-e8fee9dafefbfdc8.png?imageMogr2/auto-orient/strip|imageView2/2/w/456/format/webp)

GitHub连接

## 新建一个一个工程

新建一个工程，选择`New Directory`

![img](https:////upload-images.jianshu.io/upload_images/1716465-099c9c7de5c507d7.png?imageMogr2/auto-orient/strip|imageView2/2/w/456/format/webp)

新建工程

然后勾选`Create a git repository`

![img](https:////upload-images.jianshu.io/upload_images/1716465-907891b9916d64ef.png?imageMogr2/auto-orient/strip|imageView2/2/w/456/format/webp)

创建Git

这个时候`Rstudio`会出现`git`栏，提交到本地，只需要在`git`栏下面点击`commit`，即可提交至本地

![img](https:////upload-images.jianshu.io/upload_images/1716465-40953fcbfc448ff6.png?imageMogr2/auto-orient/strip|imageView2/2/w/456/format/webp)

GitHub提交



可以将代码保存至`GitHub`上，并且创建分支，在`GitHub`上创建一个`New respository`，命名为`test`

![img](https:////upload-images.jianshu.io/upload_images/1716465-dc1a9c7624d2109f.png?imageMogr2/auto-orient/strip|imageView2/2/w/456/format/webp)

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

![img](https:////upload-images.jianshu.io/upload_images/1716465-1e8cc2c8dfcaf26b.png?imageMogr2/auto-orient/strip|imageView2/2/w/456/format/webp)

Paste_Image.png

然后在`shell`窗口输入
 `git config remote.origin.url git@github.com:ewenharrison/test.git`

## git中设置上游

在`git`的时候，我们会建立许多有特性的分支，建立分支的时候，如何使得远端也出现分支，需要用到下面的命令：



```bash
 git push --set-upstream origin master
```



作者：cheng2pj
链接：https://www.jianshu.com/p/aa9b22b429ee
来源：简书
著作权归作者所有。商业转载请联系作者获得授权，非商业转载请注明出处。