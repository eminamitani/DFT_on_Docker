FROM ubuntu:latest
RUN apt update
RUN apt install -y build-essential
RUN apt install -y python3
RUN apt install -y gfortran
RUN apt install -y openssh-server

RUN mkdir /var/run/sshd
#SSHの設定
RUN sed -i 's/PermitRootLogin prohibit-password/PermitRootLogin yes/' /etc/ssh/sshd_config
RUN sed -i 's/#PasswordAuthentication yes/PasswordAuthentication no/' /etc/ssh/sshd_config

# SSH login fix. Otherwise user is kicked off after login
RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd

ENV NOTVISIBLE "in users profile"
RUN echo "export VISIBLE=now" >> /etc/profile

# 手元の公開鍵をコピー
COPY id_rsa_container.pub /root/authorized_keys

# ssh用の port を晒す
EXPOSE 22

# 公開鍵を使えるようにする (パーミッション変更など)
CMD mkdir ~/.ssh && \
    mv ~/authorized_keys ~/.ssh/authorized_keys && \
    chmod 0600 ~/.ssh/authorized_keys &&  \
    # 最後に ssh を起動
    /usr/sbin/sshd -D