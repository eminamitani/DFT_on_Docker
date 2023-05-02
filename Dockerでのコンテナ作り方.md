# Dockerを使って第一原理計算をお手軽に試せる仮想環境を作る

第一原理計算をためしたり、そのデータ処理をしたりするのに、一番適しているのはUbuntuに代表されるLinux環境である。しかし学部生や修士課程で自前のLinux環境を持っている人は少ないだろう。
そんなときに使えるのが仮想環境である。
WindowsやMacOS上でLinuxの仮想環境を作る方法にはいくつかあるが、
慣れてしまえばDockerが手っ取り早い。

というわけで、今回はDockerを使って仮想環境を作ってみようと思う。

## Dockerインストール

### 公式サイト
https://www.docker.com/
ここから、自分のシステムに沿ったものをダウンロードしてインストールします。
Dockerデスクトップを起動すればチュートリアルが始まるので、それが動くか確認しておきましょう。

## 公開鍵暗号
Dockerで立ち上げたコンテナ上であれこれ作業するとき、もちろん
```
docker exec -it hogehoge /bin/bash
```
のように、コンテナのbashをdockerから呼んできて作業もできますが、
ssh(公開鍵暗号)でつないで処理するほうがやりやすいです。SCPでファイル転送もできます。
というわけで、コンテナとやり取りするための公開鍵を作ってきます。
MacOSだと
```
ssh-keygen
```
でコマンド打てば秘密鍵（手元においておくほう）
と、公開鍵(.pubの拡張子がついている方)のペアが作成されます。
Windowsの場合はputty-genをインストールしてきて使う。

仮想環境なので、パスフレーズはなしで作っておく。

## DockerFileを使って必要なもの一式が入ったイメージを構築
素のUbuntuだとコンパイルとかに必要なものが全く入っていないし、
上述の公開鍵暗号通信のライブラリも入っていない。

なので、Ubuntuの最新版をベースに必要なものを予め積み込んだイメージを作っておく。
このイメージ作成方法を指定したファイルはDockerFileという。
ファイル名もDockerFileにしておく。
今回使うDockerFileの中身は以下のようなものである。

```
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
```

## イメージの構築
DockerFileを作業用ディレクトリにおいて、
-t のオプションでimageの名前を決めてビルドする。
```
docker build . -t ubuntu_ssh
```
できてるかは
```
docker image ls
```
で確認してみよう。

## コンテナの起動

```
docker run -dit --name ubuntu-2 -p 12222:22 ubuntu_ssh
```
のようにしてコンテナを起動する。`-dit`はバックグラウンドで実行＆端末からキーボード入出力を受け取ったりするためのオプション類である。
`-p 12222:22`はコンテナのport 22をホスト（自分のPC）のポート10022に連結することを指している。
このコマンドで、さっきビルドしたubuntu_sshのイメージに基づいた、コンテナubuntu-2が起動する。

作ったコンテナにssh接続してみよう。作成した秘密鍵の名前がid_rsa_containerだとすると
```
ssh root@localhost -p 12222 -i id_rsa_container
```
でssh 接続される。

いちいちコマンドを書くのが面倒な場合、MacOSであれば
`~/.ssh/config`に
```
Host docker-ubuntu
     HostName localhost
     Port 12222
     User root
     IdentityFile ~/Work/Docker/id_rsa_container
```
のように記載しておけば、
```
ssh docker-ubuntu
```
だけで接続できる。

上の手順で起動したDockerコンテナを停止するには
```
docker stop ubuntu-2
```

再び起動するには
```
docker start ubuntu-2
```



## Quantum-Espresso(q-e)のインストール
q-eはgitを使ってダウンロードしてくるのが最も早いので、上記の方法でssh接続した状態で
```
apt install git
```
を行う。

その後、q-e用のディレクトリを作り、そこでgit cloneを行う
```
mkdir q-e
cd q-e
git clone https://github.com/QEF/q-e.git .
```
とする。
これでq-eのソースコードが入手できた。
gitでダウンロードしてきたものは、現在進行系で様々な改良が加えられている(=未報告のバグがあるかも)
なので、安定バージョンになるようにgit でcheckoutする
```
git fetch
git checkout qe-7.0
```
この状態で電子状態計算の部分を担うpw.xをコンパイルしよう

```
./configure
make pw
```

## 電子状態計算を試してみる
quantum-espressoのような平面波基底のコードで最適な系というわけではないのだけれども、感覚がわかりやすいので、単純な分子の計算をしてみようと思う。

ベンゼン分子の計算に使うインプットは以下のようなものだ
```
 &control
    calculation  = 'scf'
    restart_mode = 'from_scratch'
    prefix       = 'benz'
    pseudo_dir   = './'
    outdir       = './'
 /
 &system    
    ibrav     =  1
    celldm(1) = 30
    nat       =  12
    ntyp      =  2
    ecutwfc   = 30.0
    ecutrho   = 240.0
    occupations = 'fixed' 
 /
 &electrons
    diagonalization = 'david'
    mixing_beta     =  0.2 
    conv_thr        =  1.0e-8
 /
ATOMIC_SPECIES
 C   12   C.pbe-rrkjus.UPF
 H    1   H.pbe-rrkjus.UPF
ATOMIC_POSITIONS {bohr}
C     15.628439779     15.000000000  15.00
C     14.309410418     17.284679795  15.00
C     11.671351697     17.284679795  15.00
C     10.352322336     15.000000000  15.00
C     11.671351697     12.715320205  15.00
C     14.309410418     12.715320205  15.00
H     17.675013988     15.000000000  15.00
H     15.333642386     19.057243606  15.00
H     10.647119729     19.057243606  15.00
H      8.305748128     15.000000000  15.00
H     10.647119729     10.942756393  15.00
H     15.333642386     10.942756393  15.00
K_POINTS {gamma}

```

一番最後のブロックはC6H12の分子構造のxyz座標を書いている。
その前の
```
ATOMIC_SPECIES
 C   12   C.pbe-rrkjus.UPF
 H    1   H.pbe-rrkjus.UPF
```
のところの、`C.pbe-rrkjus.UPF`は、計算に使う擬ポテンシャルの種類を指している。
http://pseudopotentials.quantum-espresso.org/legacy_tables
のWebサイトから、各原子のページに飛んで同じ名前のファイルをダウンロードし、
scpでコンテナにファイルを転送する。
SCPはFileZillaやWinSCPのようなファイルクライアントを使っても良いし、
```
scp -i id_rsa_container -P 10022 C.pbe-rrkjus.UPF root@localhost:~/q-e-test/benzen/
```
のようなコマンドで送ることもできる（コンテナ内ではq-e-test/benzenの中で作業している）

作業ディレクトリにインプットファイル(pw.inとしておく)、C.pbe-rrkjus.UPF、H.pbe-rrkjus.UPFの3つが揃ったら

```
~/q-e/bin/pw.x -in pw.in
```
とすると計算が開始される。メモリは1.5GBぐらい消費するので、計算中はPCが重くなるかもしれない。