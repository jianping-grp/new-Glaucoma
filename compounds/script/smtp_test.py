#-*-coding: utf-8-*-
# import smtplib
#
# sender = 'libaiqing11@163.com'
# receivers = 'libaiqing11@163.com'
#
# message = """From: From Person <libaiqing11@163.com>
# To: To Person <libaiqing11@163.com>
# Subject: SMTP e-mail test
#
# This is a test e-mail message.
# """
# try:
#     smtObj = smtplib.SMTP('smtp.qq.com', 587)
#     smtObj.sendmail(sender, receivers, message)
# except smtplib.SMTPException:
#     print("Error: unable to send email")

# import smtplib
# from email.mime.text import MIMEText
# from email.utils import formataddr
#
# my_sender = '1032847174@qq.com'
# my_pass = 'kyvixkkxcfadbccd'
# my_user = '1032847174@qq.com'
#
# def mail():
#     ret = True
#     try:
#         msg = MIMEText('填写邮件内容','plain', 'utf-8')
#         msg['From'] = formataddr(["FromRunoob", my_sender])  # 括号里的对应发件人邮箱昵称、发件人邮箱账号
#         msg['To'] = formataddr(["FK", my_user])              # 括号里的对应收件人邮箱昵称、收件人邮箱账号
#         msg['Subject'] = "菜鸟教程发送邮件测试"                # 邮件的主题，也可以说是标题
#
#         server = smtplib.SMTP_SSL("smtp.qq.com", 465)  # 发件人邮箱中的SMTP服务器，端口是25
#         server.login(my_sender, my_pass)  # 括号中对应的是发件人邮箱账号、邮箱密码
#         server.sendmail(my_sender, [my_user, ], msg.as_string())  # 括号中对应的是发件人邮箱账号、收件人邮箱账号、发送邮件
#         server.quit()  # 关闭连接
#     except Exception:  # 如果 try 中的语句没有执行，则会执行下面的 ret=False
#         ret = False
#     return ret
#
# ret = mail()
# if ret:
#     print("邮件发送成功")
# else:
#     print("邮件发送失败")

import smtplib
from email.mime.text import MIMEText
from email.header import Header
from smtplib import SMTP_SSL

from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.header import Header

# qq邮箱smtp服务器
host_server = 'smtp.qq.com'
# sender_qq为发件人的qq号码
sender_qq = '1032847174@qq.com'
# pwd为qq邮箱的授权码
pwd = 'kyvixkkxcfadbccd'  ##
# 发件人的邮箱
sender_qq_mail = '1032847174@qq.com'
# 收件人邮箱
receiver = '1032847174@qq.com'

# 邮件的正文内容
mail_content = "你好，<p>这是使用python登录qq邮箱发送HTML格式邮件的测试：</p><p><a href='http://www.yiibai.com'>易百教程</a></p>"
# 邮件标题
mail_title = 'Maxsu的邮件'

# 邮件正文内容
msg = MIMEMultipart()
# msg = MIMEText(mail_content, "plain", 'utf-8')
msg["Subject"] = Header(mail_title, 'utf-8')
msg["From"] = sender_qq_mail
msg["To"] = Header("接收者测试", 'utf-8')  ## 接收者的别名

# 邮件正文内容
msg.attach(MIMEText(mail_content, 'html', 'utf-8'))

# 构造附件1，传送当前目录下的 test.txt 文件
att1 = MIMEText(open('1032847174@qq.com.csv', 'rb').read(), 'base64', 'utf-8')
att1["Content-Type"] = 'application/octet-stream'
# 这里的filename可以任意写，写什么名字，邮件中显示什么名字
att1["Content-Disposition"] = 'attachment; filename="attach.txt"'
msg.attach(att1)


smtp = SMTP_SSL(host_server)
# set_debuglevel()是用来调试的。参数值为1表示开启调试模式，参数值为0关闭调试模式
smtp.set_debuglevel(1)
smtp.ehlo(host_server)
smtp.login(sender_qq, pwd)

smtp.sendmail(sender_qq_mail, receiver, msg.as_string())
smtp.quit()
