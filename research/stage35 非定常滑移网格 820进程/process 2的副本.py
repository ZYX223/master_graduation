Module global
    use mpi
    implicit none  !设计任何和隐含说明语句无效，这个时候所有变量都要显式地人工声明，不能未声明就直接使用，有效地避免了可能的大量错误
    save
    integer ::nb,nt,ng,nxm1,nym1,nzm1,nng0,nitt,nng,ign,i,j,k,it00,lbb,nxm,nym,nzm,ibm,itm,jbm,jtm
    integer,allocatable ::cg(:)
    integer ::nx,ny,nz,ib,it,jb,jt,n,ci,cycle_by
    real(8) ::vxx,vrr,vtt,vee,y1,z1,rr,sir,cor,vx,vy,vz,wx,wy,wz,ama,dim,en,pp,rpm,ma
    real(8),parameter ::pi=3.141592654
    real(8) ::c2,cfl,a4,a2,beta1,beta2
    real(8) ::time_begin,time_end
    character(len=50):: zonename,nnspan
    integer,allocatable ::nnx(:),nny(:),nnz(:),nib(:),nit(:),njb(:),njt(:)  !allocatable可变大小的数组
    real(8),allocatable ::Q01(:,:,:),Q02(:,:,:),Q03(:,:,:),Q04(:,:,:),Q05(:,:,:),Q06(:,:,:)
    real(8),allocatable ::Q31(:,:,:,:),Q32(:,:,:,:),Q33(:,:,:,:),Q34(:,:,:,:),Q35(:,:,:,:),Q36(:,:,:,:)
    real(8),allocatable ::Q11(:,:,:,:),Q12(:,:,:,:),Q13(:,:,:,:),Q14(:,:,:,:),Q15(:,:,:,:),Q16(:,:,:,:)
    real(8),allocatable ::AV1(:,:,:),AV2(:,:,:),AV3(:,:,:),AV4(:,:,:),AV5(:,:,:),AV6(:,:,:)
    real(8),allocatable ::qc1(:,:,:),qc2(:,:,:),qc3(:,:,:),qc4(:,:,:),qc5(:,:,:),qc6(:,:,:)
    real(8),allocatable ::qv1(:,:,:),qv2(:,:,:),qv3(:,:,:),qv4(:,:,:),qv5(:,:,:),qv6(:,:,:)
    real(8),allocatable ::gradfi(:,:,:,:),gradfj(:,:,:,:),gradfk(:,:,:,:),gradc(:,:,:,:),gradcs(:,:,:,:)
    real(8),allocatable ::ts1(:,:,:),ts2(:,:,:),ts3(:,:,:),ts4(:,:,:),ts5(:,:,:),ts6(:,:,:)
    real(8),allocatable ::py1(:,:,:),py2(:,:,:),py3(:,:,:),py4(:,:,:),py5(:,:,:),py6(:,:,:)
    real(8),allocatable ::rr1(:,:,:),rr2(:,:,:),rr3(:,:,:),rr4(:,:,:),rr5(:,:,:)
    real(8),allocatable ::qp1(:,:,:),qp2(:,:,:),qp3(:,:,:),qp4(:,:,:),qp5(:,:,:)
    real(8),allocatable ::xf(:,:,:),yf(:,:,:),zf(:,:,:),x(:,:,:),y(:,:,:),z(:,:,:),xx00(:,:,:),yy00(:,:,:),zz00(:,:,:)
    real(8),allocatable ::xx1(:,:,:,:),yy1(:,:,:,:),zz1(:,:,:,:),xx2(:,:,:,:),yy2(:,:,:,:),zz2(:,:,:,:),xx3(:,:,:,:),yy3(:,:,:,:),zz3(:,:,:,:),xx(:,:,:,:),yy(:,:,:,:),zz(:,:,:,:)
    real(8),allocatable ::xx01(:,:,:),yy01(:,:,:),zz01(:,:,:),xx02(:,:,:),yy02(:,:,:),zz02(:,:,:),xx03(:,:,:),yy03(:,:,:),zz03(:,:,:),xx0(:,:,:),yy0(:,:,:),zz0(:,:,:)
    real(8),allocatable ::s1xn(:,:,:,:),s1yn(:,:,:,:),s1zn(:,:,:,:),s2xn(:,:,:,:),s2yn(:,:,:,:),s2zn(:,:,:,:),s3xn(:,:,:,:),s3yn(:,:,:,:),s3zn(:,:,:,:),vvn(:,:,:,:)
    real(8),allocatable ::s1x(:,:,:),s1y(:,:,:),s1z(:,:,:),s2x(:,:,:),s2y(:,:,:),s2z(:,:,:),s3x(:,:,:),s3y(:,:,:),s3z(:,:,:),vv(:,:,:)
    real(8),allocatable ::pvx(:,:,:),pvy(:,:,:),pvz(:,:,:),vth(:,:,:),vre(:,:,:),p(:,:,:),t(:,:,:),time(:,:,:),wma(:,:,:)
    real(8),allocatable ::dmini(:,:,:)
    real(8),allocatable ::sri(:,:,:),srj(:,:,:),srk(:,:,:)
    real(8),allocatable ::dm(:,:)
    real(8),allocatable ::rms(:)
    real(8),allocatable ::betaxn(:,:,:,:),betayn(:,:,:,:),betazn(:,:,:,:),petn(:,:,:,:),pebn(:,:,:,:),turin(:,:,:,:),hatn(:,:,:,:)
    real(8),allocatable ::betax(:,:),betay(:,:),betaz(:,:),pet(:,:),peb(:,:),turi(:,:),hat(:,:)
    real(8) ::ta,timl,pt,ht,rout,pb0,pb1,period,rmsm,rmsm0,rmsmmax
    real(8) ::cvl0,t0,ts,cp,prt,prl,rg,cv1,cv2,kap,sigmav,cb1,cb2,cw1,cw2,cw3,cr1,cr2,cr3
    character(len=9):: id_m
    integer  ::myid,myidl,myidr,numprocs,ierr,rc,status(MPI_STATUS_SIZE)
    End module

    program main
    use global
    implicit none

    call MPI_INIT( ierr ) !用来初始化MPI执行环境，建立多个MPI进程之间的联系，为后续通信做准备。rc和ierr分别用来得到MPI过程调用结束后的返回结果和可能的出错信息
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr ) !标识各个MPI进程的，告诉调用该函数的进程“我是谁？”。整型变量myid和numprocs分别用来记录某一个并行执行的进程的标识和所有参加计算的进程的个数
    call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr ) !标识相应进程组中有多少个进程。得到有多少个进程参与运算，并放在numprocs中
    time_begin=MPI_Wtime() !计算逝去的时间
    call ini
    call allocation
    call shu
    call ini2
    call dodd
    call distance
    call fmg(1)
    call MPI_FINALIZE(rc) !结束MPI执行环境
    end

  subroutine ini  !读入控制参数和网格
    use global
    implicit none
    real(8) ::temp,cor1,sir1,t1,t2
    !***湍流计算所需参数
    cvl0=1.7161D-5
    t0=273.16d0
    ts=110.4d0
    rg=287.d0
    cp=rg*1.4d0/0.4d0
    prl=0.72d0
    prt=0.9d0
    cv1=7.1d0
    cv2=5.d0
    kap=0.41d0
    cb1=0.1355d0
    cb2=0.622d0
    sigmav=2.d0/3.d0
    cw1=cb1/kap**2+(1.d0+cb2)/sigmav
    cw2=0.3d0
    cw3=2.d0
    cr1=1.d0
    cr2=2.d0
    cr3=1.d0
     open(10,file="ini3.dat") !10为文件编号
     read(10,*) nt,ng,beta1,beta2,cfl,a2,a4,ht,pt,pb1,c2,rmsm0
     allocate(cg(ng)) !ng大小的数组
     do i=1,ng
        read(10,*) cg(i)
     end do
     read(10,*) nxm,nym,nzm,ibm,itm,jbm,jtm,lbb,rpm,ma
     rpm=rpm*pi/30.d0
     close(10)
     allocate(xf(nxm+1,nym+1,nzm+1))
     allocate(yf(nxm+1,nym+1,nzm+1))
     allocate(zf(nxm+1,nym+1,nzm+1))
         open(55,file="grid.dat")
         read(55,*) nxm1,nym1,nzm1
         do k=1,nzm1
            do j=1,nym1
                do i=1,nxm1
                    read(55,*) xf(i,j,k),yf(i,j,k),zf(i,j,k)!下标从1开始，先列后行存储
                end do
            end do
         end do
         close(55)

         nx=nxm
         ny=nym
         nz=nzm
         temp=2.d0/dble(lbb)*pi*myid !通道数量,将数据转换为双精度，2.d0双精度版*10的0次方
         cor1=cos(temp)
         sir1=sin(temp)
         do k=1,nz+1
             do j=1,ny+1
                 do i=1,nx+1
                     t1=yf(i,j,k)
                     t2=zf(i,j,k)
                     yf(i,j,k)=t1*cor1-t2*sir1 !切线
                     zf(i,j,k)=t1*sir1+t2*cor1 !法向
                 end do
             end do
         end do
     write(id_m,'(i3)')myid !
    end subroutine ini

subroutine allocation  !申请内存
    use global
    implicit none

    allocate(nnx(ng))
    allocate(nny(ng))
    allocate(nnz(ng))
    allocate(nib(ng))
    allocate(nit(ng))
    allocate(njb(ng))
    allocate(njt(ng))
    allocate(dm(nt,nt))
    allocate(q11(0:nx+1,0:ny+1,0:nz+1,nt))
    allocate(q12(0:nx+1,0:ny+1,0:nz+1,nt))
    allocate(q13(0:nx+1,0:ny+1,0:nz+1,nt))
    allocate(q14(0:nx+1,0:ny+1,0:nz+1,nt))
    allocate(q15(0:nx+1,0:ny+1,0:nz+1,nt))
    allocate(q16(0:nx+1,0:ny+1,0:nz+1,nt))
    allocate(q31(1:nx,1:ny,1:nz,nt))
    allocate(q32(1:nx,1:ny,1:nz,nt))
    allocate(q33(1:nx,1:ny,1:nz,nt))
    allocate(q34(1:nx,1:ny,1:nz,nt))
    allocate(q35(1:nx,1:ny,1:nz,nt))
    allocate(q36(1:nx,1:ny,1:nz,nt))
    allocate(x(1:nx+1,1:ny+1,1:nz+1))
    allocate(y(1:nx+1,1:ny+1,1:nz+1))
    allocate(z(1:nx+1,1:ny+1,1:nz+1))
    allocate(xx00(1:nx,1:ny,1:nz))
    allocate(yy00(1:nx,1:ny,1:nz))
    allocate(zz00(1:nx,1:ny,1:nz))
    allocate(dmini(nx,ny,nz))
    allocate(betaxn(ny,nz,ng,nt))
    allocate(betayn(ny,nz,ng,nt))
    allocate(betazn(ny,nz,ng,nt))
    allocate(hatn(ny,nz,ng,nt))
    allocate(petn(ny,nz,ng,nt))
    allocate(turin(ny,nz,ng,nt))
    allocate(pebn(ny,nz,ng,nt))
    allocate(betax(ny,nz))
    allocate(betay(ny,nz))
    allocate(betaz(ny,nz))
    allocate(hat(ny,nz))
    allocate(pet(ny,nz))
    allocate(turi(ny,nz))
    allocate(peb(ny,nz))
    allocate(ts1(nx,ny,nz))
    allocate(ts2(nx,ny,nz))
    allocate(ts3(nx,ny,nz))
    allocate(ts4(nx,ny,nz))
    allocate(ts5(nx,ny,nz))
    allocate(ts6(nx,ny,nz))
    allocate(s2xn(nx+1,ny,nz,ng))
    allocate(s2yn(nx+1,ny,nz,ng))
    allocate(s2zn(nx+1,ny,nz,ng))
    allocate(s3xn(nx,ny+1,nz,ng))
    allocate(s3yn(nx,ny+1,nz,ng))
    allocate(s3zn(nx,ny+1,nz,ng))
    allocate(s1xn(nx,ny,nz+1,ng))
    allocate(s1yn(nx,ny,nz+1,ng))
    allocate(s1zn(nx,ny,nz+1,ng))
    allocate(vvn(nx,ny,nz,ng))
    allocate(xx2(nx+1,ny,nz,ng))
    allocate(yy2(nx+1,ny,nz,ng))
    allocate(zz2(nx+1,ny,nz,ng))
    allocate(xx3(nx,ny+1,nz,ng))
    allocate(yy3(nx,ny+1,nz,ng))
    allocate(zz3(nx,ny+1,nz,ng))
    allocate(xx1(nx,ny,nz+1,ng))
    allocate(yy1(nx,ny,nz+1,ng))
    allocate(zz1(nx,ny,nz+1,ng))
    allocate(xx(nx,ny,nz,ng))
    allocate(yy(nx,ny,0:nz+1,ng))
    allocate(zz(nx,ny,0:nz+1,ng))
    allocate(s1x(nx,ny,0:nz+2))
    allocate(s1y(nx,ny,0:nz+2))
    allocate(s1z(nx,ny,0:nz+2))
    allocate(s2x(0:nx+2,ny,nz))
    allocate(s2y(0:nx+2,ny,nz))
    allocate(s2z(0:nx+2,ny,nz))
    allocate(s3x(nx,0:ny+2,nz))
    allocate(s3y(nx,0:ny+2,nz))
    allocate(s3z(nx,0:ny+2,nz))
    allocate(vv(nx,ny,nz))
    allocate(xx01(nx,ny,nz+1))
    allocate(yy01(nx,ny,nz+1))
    allocate(zz01(nx,ny,nz+1))
    allocate(xx02(nx+1,ny,nz))
    allocate(yy02(nx+1,ny,nz))
    allocate(zz02(nx+1,ny,nz))
    allocate(xx03(nx,ny+1,nz))
    allocate(yy03(nx,ny+1,nz))
    allocate(zz03(nx,ny+1,nz))
    allocate(xx0(0:nx+1,0:ny+1,0:nz+1))
    allocate(yy0(0:nx+1,0:ny+1,0:nz+1))
    allocate(zz0(0:nx+1,0:ny+1,0:nz+1))
    allocate(q01(1:nx,1:ny,1:nz))
    allocate(q02(1:nx,1:ny,1:nz))
    allocate(q03(1:nx,1:ny,1:nz))
    allocate(q04(1:nx,1:ny,1:nz))
    allocate(q05(1:nx,1:ny,1:nz))
    allocate(q06(1:nx,1:ny,1:nz))
    allocate(av1(0:nx+1,0:ny+1,0:nz+1))
    allocate(av2(0:nx+1,0:ny+1,0:nz+1))
    allocate(av3(0:nx+1,0:ny+1,0:nz+1))
    allocate(av4(0:nx+1,0:ny+1,0:nz+1))
    allocate(av5(0:nx+1,0:ny+1,0:nz+1))
    allocate(av6(0:nx+1,0:ny+1,0:nz+1))
    allocate(qc1(0:nx+1,0:ny+1,0:nz+1))
    allocate(qc2(0:nx+1,0:ny+1,0:nz+1))
    allocate(qc3(0:nx+1,0:ny+1,0:nz+1))
    allocate(qc4(0:nx+1,0:ny+1,0:nz+1))
    allocate(qc5(0:nx+1,0:ny+1,0:nz+1))
    allocate(qc6(0:nx+1,0:ny+1,0:nz+1))
    allocate(qv1(0:nx+1,0:ny+1,0:nz+1))
    allocate(qv2(0:nx+1,0:ny+1,0:nz+1))
    allocate(qv3(0:nx+1,0:ny+1,0:nz+1))
    allocate(qv4(0:nx+1,0:ny+1,0:nz+1))
    allocate(qv5(0:nx+1,0:ny+1,0:nz+1))
    allocate(qv6(0:nx+1,0:ny+1,0:nz+1))
    allocate(gradfi(15,1:nx+1,0:ny+1,0:nz+1))
    allocate(gradfj(15,0:nx+1,1:ny+1,0:nz+1))
    allocate(gradfk(15,0:nx+1,0:ny+1,1:nz+1))
    allocate(gradc(12,0:nx+1,0:ny+1,0:nz+1))
    allocate(gradcs(9,0:nx+1,0:ny+1,0:nz+1))
    allocate(py1(nx,ny,nz))
    allocate(py2(nx,ny,nz))
    allocate(py3(nx,ny,nz))
    allocate(py4(nx,ny,nz))
    allocate(py5(nx,ny,nz))
    allocate(py6(nx,ny,nz))
    allocate(rr1(nx,ny,nz))
    allocate(rr2(nx,ny,nz))
    allocate(rr3(nx,ny,nz))
    allocate(rr4(nx,ny,nz))
    allocate(rr5(nx,ny,nz))
    allocate(qp1(nx,ny,nz))
    allocate(qp2(nx,ny,nz))
    allocate(qp3(nx,ny,nz))
    allocate(qp4(nx,ny,nz))
    allocate(qp5(nx,ny,nz))
    allocate(pvx(0:nx+1,0:ny+1,0:nz+1))
    allocate(pvy(0:nx+1,0:ny+1,0:nz+1))
    allocate(pvz(0:nx+1,0:ny+1,0:nz+1))
    allocate(vth(0:nx+1,0:ny+1,0:nz+1))
    allocate(vre(0:nx+1,0:ny+1,0:nz+1))
    allocate(p(0:nx+1,0:ny+1,0:nz+1))
    allocate(t(0:nx+1,0:ny+1,0:nz+1))
    allocate(wma(0:nx+1,0:ny+1,0:nz+1))
    allocate(time(nx,ny,nz))
    allocate(rms(nt))
    allocate(sri(0:nx+1,0:ny+1,0:nz+1))
    allocate(srj(0:nx+1,0:ny+1,0:nz+1))
    allocate(srk(0:nx+1,0:ny+1,0:nz+1))
    end subroutine allocation

  subroutine shu  !多重网格数
    use global
    implicit none
    integer :: kk

        itm=itm-1 !叶片处的网格，itm已读入
        jtm=jtm-1
        do j=1,ng
            kk=2**(ng-j) !网格数量，粗，中间，细
            nnx(j)=nx/kk
            nny(j)=ny/kk
            nnz(j)=nz/kk
            nib(j)=(ibm-1)/kk+1 !叶片处的网格
            nit(j)=itm/kk       !
            njb(j)=(jbm-1)/kk+1
            njt(j)=jtm/kk
        end do
    end subroutine shu

subroutine ini2 !最粗网格上初始化流场
    use global
    implicit none
    real(8) ::cvl,a
    integer ::kk,ifine,jfine,kfine

    nng0=1 !当前网格层，其中粗网格层为1，中网格层为2，细为3
    nx=nnx(nng0)
    ny=nny(nng0)
    nz=nnz(nng0)
    kk=2**(ng-nng0)!细化率
    do k=1,nz+1
        kfine=kk*(k-1)+1
        do j=1,ny+1
            jfine=kk*(j-1)+1
            do i=1,nx+1
                ifine=kk*(i-1)+1
                x(i,j,k)=xf(ifine,jfine,kfine) !粗网格上顶点处坐标，隔4取1个
                y(i,j,k)=yf(ifine,jfine,kfine)
                z(i,j,k)=zf(ifine,jfine,kfine)
            end do
        end do
    end do
    vxx=1.d0 !求气流角时用到的，三个方向
    vtt=vxx*beta1
    vrr=vxx*beta2
    vee=sqrt(vxx*vxx+vrr*vrr+vtt*vtt)
    do k=1,nz
        do j=1,ny
            y1=0.25d0*(y(1,j,k)+y(1,j+1,k)+y(1,j,k+1)+y(1,j+1,k+1)) !x面上的值
            z1=0.25d0*(z(1,j,k)+z(1,j+1,k)+z(1,j,k+1)+z(1,j+1,k+1))
            rr=sqrt(y1*y1+z1*z1)
            sir=z1/rr
            cor=y1/rr
            betaxn(j,k,nng0,1:nt)=vxx/vee
            betayn(j,k,nng0,1:nt)=(vrr*cor-vtt*sir)/vee
            betazn(j,k,nng0,1:nt)=(vrr*sir+vtt*cor)/vee
        end do
    end do
    hatn(1:ny,1:nz,nng0,1:nt)=ht  !进口总焓，ht为之前读文件得到的值
    petn(1:ny,1:nz,nng0,1:nt)=pt  !进口总压，pt为之前读文件得到的值
        nitt=0   !总迭代次数
        do n=1,nt
            do k=1,nz
                do j=1,ny
                    do i=1,nx
                        rout=3.5d0*petn(j,k,nng0,n)/hatn(j,k,nng0,n)
                        pp=petn(j,k,nng0,n)*(1.d0+0.2d0*ma*ma)**(-3.5)
                        dim=rout*(1.d0+0.2d0*ma*ma)**(-2.5)
                        a=sqrt(1.4d0*pp/dim)
                        en=2.5d0*pp+0.5d0*dim*ma*ma*a*a
                        q11(i,j,k,n)=dim !原始守衡量
                        q12(i,j,k,n)=dim*a*ma*betaxn(j,k,nng0,n)
                        q13(i,j,k,n)=dim*a*ma*betayn(j,k,nng0,n)
                        q14(i,j,k,n)=dim*a*ma*betazn(j,k,nng0,n)
                        q15(i,j,k,n)=en
                        q16(i,j,k,n)=200.d0*cvl0
                    end do
                end do
            end do
        end do
        open(myid+200,file='error-'//trim(adjustl(id_m))//'myid.dat') !打开文件
        if(myid==0)then
            open(404+myid,file='convergence-inflow.dat')
            open(707+myid,file='convergence-outflow.dat')
        end if
    end subroutine ini2

subroutine dodd!时间谱系数
    use global
    implicit none

    period=2.d0*pi/rpm
    dm=0.d0
    if(mod(nt,2)==1) then
        do j=1,nt  !nt=1?
            do i=1,nt
                if (i/=j) then
                    dm(i,j)=pi/period*(-1.d0)**(i-j)/sin(pi*dble((i-j))/dble(nt)) !双精度实型二维数组；时间谱的加权因子
                end if
            end do
        end do
    else
        do j=1,nt
            do i=1,nt
                if (i/=j) then
                    dm(i,j)=pi/period*(-1.d0)**(i-j)/tan(pi*dble((i-j))/dble(nt))
                end if
            end do
        end do
    end if
    end subroutine dodd

subroutine distance    !计算最细网格上到壁面的最短距离
    use global
    implicit none
    integer ::xb,xt,ii,jj,kk,iil,iir,kkl,kkr,iib,iit,jjb,jjt
    real(8) ::d,dy1,dy2,dy,dz1,dz2,dz
    real(8),allocatable ::xwu(:,:),ywu(:,:),zwu(:,:),xwd(:,:),ywd(:,:),zwd(:,:),xwf(:,:),ywf(:,:),zwf(:,:),xwb(:,:),ywb(:,:),zwb(:,:)
    nx=nnx(ng)
    ny=nny(ng)
    nz=nnz(ng)
    ib=nib(ng)
    it=nit(ng)
    jb=njb(ng)
    jt=njt(ng)
        do k=1,nz
            do j=1,ny
                do i=1,nx !中心点处坐标
                    xx00(i,j,k)=0.125d0*(xf(i,j,k)+xf(i+1,j,k)+xf(i,j,k+1)+xf(i+1,j,k+1)+xf(i,j+1,k)+xf(i+1,j+1,k)+xf(i,j+1,k+1)+xf(i+1,j+1,k+1))
                    yy00(i,j,k)=0.125d0*(yf(i,j,k)+yf(i+1,j,k)+yf(i,j,k+1)+yf(i+1,j,k+1)+yf(i,j+1,k)+yf(i+1,j+1,k)+yf(i,j+1,k+1)+yf(i+1,j+1,k+1))
                    zz00(i,j,k)=0.125d0*(zf(i,j,k)+zf(i+1,j,k)+zf(i,j,k+1)+zf(i+1,j,k+1)+zf(i,j+1,k)+zf(i+1,j+1,k)+zf(i,j+1,k+1)+zf(i+1,j+1,k+1))
                end do
            end do
        end do
        allocate(xwd(1:nx,1:nz))
        allocate(ywd(1:nx,1:nz))
        allocate(zwd(1:nx,1:nz))
        allocate(xwu(1:nx,1:nz))
        allocate(ywu(1:nx,1:nz))
        allocate(zwu(1:nx,1:nz))
        do k=1,nz  !y边界
            do i=1,nx
                xwd(i,k)=0.25d0*(xf(i,1,k)+xf(i+1,1,k)+xf(i,1,k+1)+xf(i+1,1,k+1))
                ywd(i,k)=0.25d0*(yf(i,1,k)+yf(i+1,1,k)+yf(i,1,k+1)+yf(i+1,1,k+1))
                zwd(i,k)=0.25d0*(zf(i,1,k)+zf(i+1,1,k)+zf(i,1,k+1)+zf(i+1,1,k+1))

                xwu(i,k)=0.25d0*(xf(i,ny+1,k)+xf(i+1,ny+1,k)+xf(i,ny+1,k+1)+xf(i+1,ny+1,k+1))
                ywu(i,k)=0.25d0*(yf(i,ny+1,k)+yf(i+1,ny+1,k)+yf(i,ny+1,k+1)+yf(i+1,ny+1,k+1))
                zwu(i,k)=0.25d0*(zf(i,ny+1,k)+zf(i+1,ny+1,k)+zf(i,ny+1,k+1)+zf(i+1,ny+1,k+1))
            end do
        end do
    allocate(xwf(ib:it,jb:jt))
    allocate(ywf(ib:it,jb:jt))
    allocate(zwf(ib:it,jb:jt))
    allocate(xwb(ib:it,jb:jt))
    allocate(ywb(ib:it,jb:jt))
    allocate(zwb(ib:it,jb:jt))
        do j=jb,jt
            do i=ib,it   !z
                 xwf(i,j)=0.25d0*(xf(i,j,1)+xf(i+1,j,1)+xf(i,j+1,1)+xf(i+1,j+1,1))
                 ywf(i,j)=0.25d0*(yf(i,j,1)+yf(i+1,j,1)+yf(i,j+1,1)+yf(i+1,j+1,1))
                 zwf(i,j)=0.25d0*(zf(i,j,1)+zf(i+1,j,1)+zf(i,j+1,1)+zf(i+1,j+1,1))

                 xwb(i,j)=0.25d0*(xf(i,j,nz+1)+xf(i+1,j,nz+1)+xf(i,j+1,nz+1)+xf(i+1,j+1,nz+1))
                 ywb(i,j)=0.25d0*(yf(i,j,nz+1)+yf(i+1,j,nz+1)+yf(i,j+1,nz+1)+yf(i+1,j+1,nz+1))
                 zwb(i,j)=0.25d0*(zf(i,j,nz+1)+zf(i+1,j,nz+1)+zf(i,j+1,nz+1)+zf(i+1,j+1,nz+1))
            end do
        end do

        do k=1,nz
            do j=1,ny
                do i=1,nx
                    dy=1.d20
                    iil=i-10
                    iil=max(iil,1)
                    iir=iil+20
                    iir=min(iir,nx)

                    kkl=k-10
                    kkl=max(kkl,1)
                    kkr=kkl+20
                    kkr=min(kkr,nz)

                    do kk=kkl,kkr
                        do ii=iil,iir
                            dy1=(xx00(i,j,k)-xwd(ii,kk))**2+(yy00(i,j,k)-ywd(ii,kk))**2+(zz00(i,j,k)-zwd(ii,kk))**2
                            dy2=(xx00(i,j,k)-xwu(ii,kk))**2+(yy00(i,j,k)-ywu(ii,kk))**2+(zz00(i,j,k)-zwu(ii,kk))**2
                            d=min(dy1,dy2)
                            if(d<dy)then
                                dy=d
                            end if
                        end do
                    end do

                    dz=1.d20
                    if (i<ib) then
                        iib=ib
                        iit=ib+20
                    else if(i>=ib .and. i<=it) then
                        iib=i-10
                        iib=max(iib,ib)
                        iit=iib+20
                        iit=min(iit,it)
                    else
                        iib=it-20
                        iit=it
                    end if

                    if (j<jb) then
                        jjb=jb
                        jjt=jb+20
                    else if (j>=jb .and. j<=jt) then
                        jjb=j-10
                        jjb=max(jjb,jb)
                        jjt=jjb+20
                        jjt=min(jjt,jt)
                    else
                        jjb=jt-20
                        jjt=jt
                    end if

                    do jj=jjb,jjt
                        do ii=iib,iit
                            dz1=(xx00(i,j,k)-xwf(ii,jj))**2+(yy00(i,j,k)-ywf(ii,jj))**2+(zz00(i,j,k)-zwf(ii,jj))**2
                            dz2=(xx00(i,j,k)-xwb(ii,jj))**2+(yy00(i,j,k)-ywb(ii,jj))**2+(zz00(i,j,k)-zwb(ii,jj))**2
                            d=min(dz1,dz2)
                            if (d<dz) then
                                dz=d
                            end if
                        end do
                    end do
                    dmini(i,j,k)=sqrt(min(dy,dz))  !求出最短距离
                end do
            end do
        end do
        deallocate(xwd)
        deallocate(ywd)
        deallocate(zwd)
        deallocate(xwu)
        deallocate(ywu)
        deallocate(zwu)
        deallocate(xwf)
        deallocate(ywf)
        deallocate(zwf)
        deallocate(xwb)
        deallocate(ywb)
        deallocate(zwb)
    end subroutine distance

subroutine fmg(style)
    use global
    implicit none
    integer,intent(in) :: style !intent(in)表示这个参数是输入的

    FMGcycle: do nng=nng0,ng !ng=3，遍历网格层,nng0=1
       call geo !计算网格参数
       call ye  !边界值
       Vcycle: do it00=1,cg(nng) !遍历设置好的每层网格的迭代步数
           nitt=nitt+1 !总迭代次数,此时应该是1
           rmsm=-11.d0
           Timelevel: do n=1,nt  !nt=1，遍历物理时间步数
              !*******************细网格计算****
                    call geon(nng) !具体到某层，当前网格层计算网格参数
                    call yen(nng)  !当前网格层边界条件设置
                    qp1=0.d0 !驱动源项P2h，多重网格粗网格迭代中
                    qp2=0.d0
                    qp3=0.d0
                    qp4=0.d0
                    qp5=0.d0
                    call tsd      !时间谱中的 时间导数项
                    call march1   !3步RK推进
                    call residual !残差
                    call ppp      !由守恒量计算原始量
                    call bc       !边界条件（6个面，有并行，同向数据传递）
                    call step     !当地时间步长
                    call ddd      !人工粘性
                    call qqq      !对流通量
                    call qqqv     !粘性通量
                    q01(1:nx,1:ny,1:nz)=q11(1:nx,1:ny,1:nz,n) !q01 :R-K中存放上一次时间推进结果（推进初值）, q11:原始守恒量
                    q02(1:nx,1:ny,1:nz)=q12(1:nx,1:ny,1:nz,n)
                    q03(1:nx,1:ny,1:nz)=q13(1:nx,1:ny,1:nz,n)
                    q04(1:nx,1:ny,1:nz)=q14(1:nx,1:ny,1:nz,n)
                    q05(1:nx,1:ny,1:nz)=q15(1:nx,1:ny,1:nz,n)
                    timl=cfl*0.125d0
                    call pred(0)
                    if (nng==ng) then
                        call marchsa   !****湍流模型计算***（只在最细网格上进行）
                    end if
                    q31(1:nx,1:ny,1:nz,:)=q11(1:nx,1:ny,1:nz,:) !q31:粗网格将残差插值到细网格时的修正值
                    q32(1:nx,1:ny,1:nz,:)=q12(1:nx,1:ny,1:nz,:)
                    q33(1:nx,1:ny,1:nz,:)=q13(1:nx,1:ny,1:nz,:)
                    q34(1:nx,1:ny,1:nz,:)=q14(1:nx,1:ny,1:nz,:)
                    q35(1:nx,1:ny,1:nz,:)=q15(1:nx,1:ny,1:nz,:)
                    q36(1:nx,1:ny,1:nz,:)=q16(1:nx,1:ny,1:nz,:)
                    if(style==1) then
                !********************粗网格计算*********************
                        Coarse: do ign=nng-1,1,-1
                            call geon(ign)
                            call init   !限制算子
                            call yen(ign) !
                            call copr
                            call march2
                            call update  !粗网格的残差来修正细网格，对粗网格本身无影响
                            call ppp
                            call bc
                            call step
                            call dddc  !粗网格人工粘性
                            call qqq
                            call qqqv
                            q01(1:nx,1:ny,1:nz)=q11(1:nx,1:ny,1:nz,n)
                            q02(1:nx,1:ny,1:nz)=q12(1:nx,1:ny,1:nz,n)
                            q03(1:nx,1:ny,1:nz)=q13(1:nx,1:ny,1:nz,n)
                            q04(1:nx,1:ny,1:nz)=q14(1:nx,1:ny,1:nz,n)
                            q05(1:nx,1:ny,1:nz)=q15(1:nx,1:ny,1:nz,n)
                            timl=cfl*0.1d0
                            call pred(0)
                            call update
                        end do Coarse
                        !*****************粗网格计算完成******************
                        call geon(nng)
                        q11(1:nx,1:ny,1:nz,:)=q31(1:nx,1:ny,1:nz,:)
                        q12(1:nx,1:ny,1:nz,:)=q32(1:nx,1:ny,1:nz,:)
                        q13(1:nx,1:ny,1:nz,:)=q33(1:nx,1:ny,1:nz,:)
                        q14(1:nx,1:ny,1:nz,:)=q34(1:nx,1:ny,1:nz,:)
                        q15(1:nx,1:ny,1:nz,:)=q35(1:nx,1:ny,1:nz,:)
                        q16(1:nx,1:ny,1:nz,:)=q36(1:nx,1:ny,1:nz,:)
                    end if
            end do Timelevel
            time_end=MPI_Wtime()
            write(myid+200,"(i5,2x,F10.5,x,F15.5)") nitt,rmsm,time_end-time_begin
            call MPI_ALLREDUCE(rmsm,rmsmmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
            if(rmsmmax==rmsm0) then
                stop
            end if
            call flow
        end do Vcycle
        if(nng<ng) then
            call update2
        else
            call output
        end if
    end do FMGcycle
    end subroutine fmg

subroutine geo     !计算网格参数
    use global
    implicit none
    integer ::ifine,jfine,kfine,kk
    real(8) ::val1,val2,val3,val4,val5,val6,val7,val8,val,temp,cor1,sir1
     !当前网格层
        nx=nnx(nng)
        ny=nny(nng)
        nz=nnz(nng)
        kk=2**(ng-nng)
        do k=1,nz+1
            kfine=kk*(k-1)+1
            do j=1,ny+1
                jfine=kk*(j-1)+1
                do i=1,nx+1
                    ifine=kk*(i-1)+1
                    x(i,j,k)=xf(ifine,jfine,kfine)
                    y(i,j,k)=yf(ifine,jfine,kfine)
                    z(i,j,k)=zf(ifine,jfine,kfine)
                end do
            end do
        end do
        call geo1(x(1:nx+1,1:ny+1,1:nz+1),y(1:nx+1,1:ny+1,1:nz+1),z(1:nx+1,1:ny+1,1:nz+1)  &
                 ,s1xn(1:nx,1:ny,1:nz+1,nng),s1yn(1:nx,1:ny,1:nz+1,nng),s1zn(1:nx,1:ny,1:nz+1,nng)  &
                 ,xx1(1:nx,1:ny,1:nz+1,nng),yy1(1:nx,1:ny,1:nz+1,nng),zz1(1:nx,1:ny,1:nz+1,nng)     &
                 ,s2xn(1:nx+1,1:ny,1:nz,nng),s2yn(1:nx+1,1:ny,1:nz,nng),s2zn(1:nx+1,1:ny,1:nz,nng)  &
                 ,xx2(1:nx+1,1:ny,1:nz,nng),yy2(1:nx+1,1:ny,1:nz,nng),zz2(1:nx+1,1:ny,1:nz,nng)     &
                 ,s3xn(1:nx,1:ny+1,1:nz,nng),s3yn(1:nx,1:ny+1,1:nz,nng),s3zn(1:nx,1:ny+1,1:nz,nng)  &
                 ,xx3(1:nx,1:ny+1,1:nz,nng),yy3(1:nx,1:ny+1,1:nz,nng),zz3(1:nx,1:ny+1,1:nz,nng)     &
                 ,vvn(1:nx,1:ny,1:nz,nng),xx(1:nx,1:ny,1:nz,nng),yy(1:nx,1:ny,1:nz,nng),zz(1:nx,1:ny,1:nz,nng))

        do ign=nng-1,1,-1 !每次减1，终值为1   细网格算粗网格的数据后推算前一层
            nx=nnx(ign)
            ny=nny(ign)
            nz=nnz(ign)
            do k=1,nz
                kfine=2*k-1
                do j=1,ny
                    jfine=2*j-1
                    do i=1,nx+1
                        ifine=2*i-1
                        val1=sqrt(s2xn(ifine,jfine,kfine,ign+1)**2+s2yn(ifine,jfine,kfine,ign+1)**2+s2zn(ifine,jfine,kfine,ign+1)**2) !距离
                        val2=sqrt(s2xn(ifine,jfine+1,kfine,ign+1)**2+s2yn(ifine,jfine+1,kfine,ign+1)**2+s2zn(ifine,jfine+1,kfine,ign+1)**2)
                        val3=sqrt(s2xn(ifine,jfine,kfine+1,ign+1)**2+s2yn(ifine,jfine,kfine+1,ign+1)**2+s2zn(ifine,jfine,kfine+1,ign+1)**2)
                        val4=sqrt(s2xn(ifine,jfine+1,kfine+1,ign+1)**2+s2yn(ifine,jfine+1,kfine+1,ign+1)**2+s2zn(ifine,jfine+1,kfine+1,ign+1)**2)
                        val=val1+val2+val3+val4 !xx2:面中心值，1为周向，2为轴向，3为展向
                        xx2(i,j,k,ign)=(xx2(ifine,jfine,kfine,ign+1)*val1+xx2(ifine,jfine+1,kfine,ign+1)*val2+xx2(ifine,jfine,kfine+1,ign+1)*val3+xx2(ifine,jfine+1,kfine+1,ign+1)*val4)/val
                        yy2(i,j,k,ign)=(yy2(ifine,jfine,kfine,ign+1)*val1+yy2(ifine,jfine+1,kfine,ign+1)*val2+yy2(ifine,jfine,kfine+1,ign+1)*val3+yy2(ifine,jfine+1,kfine+1,ign+1)*val4)/val
                        zz2(i,j,k,ign)=(zz2(ifine,jfine,kfine,ign+1)*val1+zz2(ifine,jfine+1,kfine,ign+1)*val2+zz2(ifine,jfine,kfine+1,ign+1)*val3+zz2(ifine,jfine+1,kfine+1,ign+1)*val4)/val

                        s2xn(i,j,k,ign)=s2xn(ifine,jfine,kfine,ign+1)+s2xn(ifine,jfine+1,kfine,ign+1)+s2xn(ifine,jfine,kfine+1,ign+1)+s2xn(ifine,jfine+1,kfine+1,ign+1)
                        s2yn(i,j,k,ign)=s2yn(ifine,jfine,kfine,ign+1)+s2yn(ifine,jfine+1,kfine,ign+1)+s2yn(ifine,jfine,kfine+1,ign+1)+s2yn(ifine,jfine+1,kfine+1,ign+1)
                        s2zn(i,j,k,ign)=s2zn(ifine,jfine,kfine,ign+1)+s2zn(ifine,jfine+1,kfine,ign+1)+s2zn(ifine,jfine,kfine+1,ign+1)+s2zn(ifine,jfine+1,kfine+1,ign+1)
                    end do
                end do
            end do

            do k=1,nz
                kfine=2*k-1
                do j=1,ny+1
                    jfine=2*j-1
                    do i=1,nx
                        ifine=2*i-1
                        val1=sqrt(s3xn(ifine,jfine,kfine,ign+1)**2+s3yn(ifine,jfine,kfine,ign+1)**2+s3zn(ifine,jfine,kfine,ign+1)**2)
                        val2=sqrt(s3xn(ifine+1,jfine,kfine,ign+1)**2+s3yn(ifine+1,jfine,kfine,ign+1)**2+s3zn(ifine+1,jfine,kfine,ign+1)**2)
                        val3=sqrt(s3xn(ifine,jfine,kfine+1,ign+1)**2+s3yn(ifine,jfine,kfine+1,ign+1)**2+s3zn(ifine,jfine,kfine+1,ign+1)**2)
                        val4=sqrt(s3xn(ifine+1,jfine,kfine+1,ign+1)**2+s3yn(ifine+1,jfine,kfine+1,ign+1)**2+s3zn(ifine+1,jfine,kfine+1,ign+1)**2)
                        val=val1+val2+val3+val4
                        xx3(i,j,k,ign)=(xx3(ifine,jfine,kfine,ign+1)*val1+xx3(ifine+1,jfine,kfine,ign+1)*val2+xx3(ifine,jfine,kfine+1,ign+1)*val3+xx3(ifine+1,jfine,kfine+1,ign+1)*val4)/val
                        yy3(i,j,k,ign)=(yy3(ifine,jfine,kfine,ign+1)*val1+yy3(ifine+1,jfine,kfine,ign+1)*val2+yy3(ifine,jfine,kfine+1,ign+1)*val3+yy3(ifine+1,jfine,kfine+1,ign+1)*val4)/val
                        zz3(i,j,k,ign)=(zz3(ifine,jfine,kfine,ign+1)*val1+zz3(ifine+1,jfine,kfine,ign+1)*val2+zz3(ifine,jfine,kfine+1,ign+1)*val3+zz3(ifine+1,jfine,kfine+1,ign+1)*val4)/val

                        s3xn(i,j,k,ign)=s3xn(ifine,jfine,kfine,ign+1)+s3xn(ifine+1,jfine,kfine,ign+1)+s3xn(ifine,jfine,kfine+1,ign+1)+s3xn(ifine+1,jfine,kfine+1,ign+1)
                        s3yn(i,j,k,ign)=s3yn(ifine,jfine,kfine,ign+1)+s3yn(ifine+1,jfine,kfine,ign+1)+s3yn(ifine,jfine,kfine+1,ign+1)+s3yn(ifine+1,jfine,kfine+1,ign+1)
                        s3zn(i,j,k,ign)=s3zn(ifine,jfine,kfine,ign+1)+s3zn(ifine+1,jfine,kfine,ign+1)+s3zn(ifine,jfine,kfine+1,ign+1)+s3zn(ifine+1,jfine,kfine+1,ign+1)
                    end do
                end do
            end do

            do k=1,nz+1
                kfine=2*k-1
                do j=1,ny
                    jfine=2*j-1
                    do i=1,nx
                        ifine=2*i-1
                        val1=sqrt(s1xn(ifine,jfine,kfine,ign+1)**2+s1yn(ifine,jfine,kfine,ign+1)**2+s1zn(ifine,jfine,kfine,ign+1)**2)
                        val2=sqrt(s1xn(ifine+1,jfine,kfine,ign+1)**2+s1yn(ifine+1,jfine,kfine,ign+1)**2+s1zn(ifine+1,jfine,kfine,ign+1)**2)
                        val3=sqrt(s1xn(ifine,jfine+1,kfine,ign+1)**2+s1yn(ifine,jfine+1,kfine,ign+1)**2+s1zn(ifine,jfine+1,kfine,ign+1)**2)
                        val4=sqrt(s1xn(ifine+1,jfine+1,kfine,ign+1)**2+s1yn(ifine+1,jfine+1,kfine,ign+1)**2+s1zn(ifine+1,jfine+1,kfine,ign+1)**2)
                        val=val1+val2+val3+val4
                        xx1(i,j,k,ign)=(xx1(ifine,jfine,kfine,ign+1)*val1+xx1(ifine+1,jfine,kfine,ign+1)*val2+xx1(ifine,jfine+1,kfine,ign+1)*val3+xx1(ifine+1,jfine+1,kfine,ign+1)*val4)/val
                        yy1(i,j,k,ign)=(yy1(ifine,jfine,kfine,ign+1)*val1+yy1(ifine+1,jfine,kfine,ign+1)*val2+yy1(ifine,jfine+1,kfine,ign+1)*val3+yy1(ifine+1,jfine+1,kfine,ign+1)*val4)/val
                        zz1(i,j,k,ign)=(zz1(ifine,jfine,kfine,ign+1)*val1+zz1(ifine+1,jfine,kfine,ign+1)*val2+zz1(ifine,jfine+1,kfine,ign+1)*val3+zz1(ifine+1,jfine+1,kfine,ign+1)*val4)/val

                        s1xn(i,j,k,ign)=s1xn(ifine,jfine,kfine,ign+1)+s1xn(ifine+1,jfine,kfine,ign+1)+s1xn(ifine,jfine+1,kfine,ign+1)+s1xn(ifine+1,jfine+1,kfine,ign+1)
                        s1yn(i,j,k,ign)=s1yn(ifine,jfine,kfine,ign+1)+s1yn(ifine+1,jfine,kfine,ign+1)+s1yn(ifine,jfine+1,kfine,ign+1)+s1yn(ifine+1,jfine+1,kfine,ign+1)
                        s1zn(i,j,k,ign)=s1zn(ifine,jfine,kfine,ign+1)+s1zn(ifine+1,jfine,kfine,ign+1)+s1zn(ifine,jfine+1,kfine,ign+1)+s1zn(ifine+1,jfine+1,kfine,ign+1)
                    end do
                end do
            end do

            do k=1,nz
                kfine=2*k-1
                do j=1,ny
                    jfine=2*j-1
                    do i=1,nx
                        ifine=2*i-1  !vvn之前并没有计算啊？
                        val1=vvn(ifine,jfine,kfine,ign+1)
                        val2=vvn(ifine+1,jfine,kfine,ign+1)
                        val3=vvn(ifine,jfine+1,kfine,ign+1)
                        val4=vvn(ifine+1,jfine+1,kfine,ign+1)
                        val5=vvn(ifine,jfine,kfine+1,ign+1)
                        val6=vvn(ifine+1,jfine,kfine+1,ign+1)
                        val7=vvn(ifine,jfine+1,kfine+1,ign+1)
                        val8=vvn(ifine+1,jfine+1,kfine+1,ign+1)
                        val=val1+val2+val3+val4+val5+val6+val7+val8
                        xx(i,j,k,ign)=(xx(ifine,jfine,kfine,ign+1)*val1+xx(ifine+1,jfine,kfine,ign+1)*val2+xx(ifine,jfine+1,kfine,ign+1)*val3+xx(ifine+1,jfine+1,kfine,ign+1)*val4+            &
                                      xx(ifine,jfine,kfine+1,ign+1)*val5+xx(ifine+1,jfine,kfine+1,ign+1)*val6+xx(ifine,jfine+1,kfine+1,ign+1)*val7+xx(ifine+1,jfine+1,kfine+1,ign+1)*val8 )/val
                        yy(i,j,k,ign)=(yy(ifine,jfine,kfine,ign+1)*val1+yy(ifine+1,jfine,kfine,ign+1)*val2+yy(ifine,jfine+1,kfine,ign+1)*val3+yy(ifine+1,jfine+1,kfine,ign+1)*val4+            &
                                      yy(ifine,jfine,kfine+1,ign+1)*val5+yy(ifine+1,jfine,kfine+1,ign+1)*val6+yy(ifine,jfine+1,kfine+1,ign+1)*val7+yy(ifine+1,jfine+1,kfine+1,ign+1)*val8 )/val
                        zz(i,j,k,ign)=(zz(ifine,jfine,kfine,ign+1)*val1+zz(ifine+1,jfine,kfine,ign+1)*val2+zz(ifine,jfine+1,kfine,ign+1)*val3+zz(ifine+1,jfine+1,kfine,ign+1)*val4+            &
                                      zz(ifine,jfine,kfine+1,ign+1)*val5+zz(ifine+1,jfine,kfine+1,ign+1)*val6+zz(ifine,jfine+1,kfine+1,ign+1)*val7+zz(ifine+1,jfine+1,kfine+1,ign+1)*val8 )/val
                        vvn(i,j,k,ign)=val
                    end do
                end do
            end do
        end do
        temp=2.d0/dble(lbb)*pi !角度，lbb周期？
        cor1=cos(temp)
        sir1=sin(temp)
            do ign=1,nng
                nx=nnx(ign)
                ny=nny(ign)
                nz=nnz(ign)
                do j=1,ny
                    do i=1,nx !z的两个虚拟边界
                        yy(i,j,nz+1,ign)=yy(i,j,1,ign)*cor1-zz(i,j,1,ign)*sir1
                        zz(i,j,nz+1,ign)=zz(i,j,1,ign)*cor1+yy(i,j,1,ign)*sir1
                        yy(i,j,0,ign)=yy(i,j,nz,ign)*cor1+zz(i,j,nz,ign)*sir1
                        zz(i,j,0,ign)=zz(i,j,nz,ign)*cor1-yy(i,j,nz,ign)*sir1
                    end do
                end do
            end do
    end subroutine geo

    subroutine geon(nnng) !具体到某层，当前网格层计算网格参数
    use global
    implicit none
    integer,intent(in) ::nnng

        nx=nnx(nnng)
        ny=nny(nnng)
        nz=nnz(nnng)
        ib=nib(nnng)
        it=nit(nnng)
        jb=njb(nnng)
        jt=njt(nnng)

        xx01(1:nx,1:ny,1:nz+1)=xx1(1:nx,1:ny,1:nz+1,nnng)
        yy01(1:nx,1:ny,1:nz+1)=yy1(1:nx,1:ny,1:nz+1,nnng)
        zz01(1:nx,1:ny,1:nz+1)=zz1(1:nx,1:ny,1:nz+1,nnng)
        xx02(1:nx+1,1:ny,1:nz)=xx2(1:nx+1,1:ny,1:nz,nnng)
        yy02(1:nx+1,1:ny,1:nz)=yy2(1:nx+1,1:ny,1:nz,nnng)
        zz02(1:nx+1,1:ny,1:nz)=zz2(1:nx+1,1:ny,1:nz,nnng)
        xx03(1:nx,1:ny+1,1:nz)=xx3(1:nx,1:ny+1,1:nz,nnng)
        yy03(1:nx,1:ny+1,1:nz)=yy3(1:nx,1:ny+1,1:nz,nnng)
        zz03(1:nx,1:ny+1,1:nz)=zz3(1:nx,1:ny+1,1:nz,nnng)
        xx0(1:nx,1:ny,1:nz)=xx(1:nx,1:ny,1:nz,nnng)
        yy0(1:nx,1:ny,0:nz+1)=yy(1:nx,1:ny,0:nz+1,nnng)
        zz0(1:nx,1:ny,0:nz+1)=zz(1:nx,1:ny,0:nz+1,nnng)

        s1x(1:nx,1:ny,1:nz+1)=s1xn(1:nx,1:ny,1:nz+1,nnng)
        s1y(1:nx,1:ny,1:nz+1)=s1yn(1:nx,1:ny,1:nz+1,nnng)
        s1z(1:nx,1:ny,1:nz+1)=s1zn(1:nx,1:ny,1:nz+1,nnng)
        s2x(1:nx+1,1:ny,1:nz)=s2xn(1:nx+1,1:ny,1:nz,nnng)
        s2y(1:nx+1,1:ny,1:nz)=s2yn(1:nx+1,1:ny,1:nz,nnng)
        s2z(1:nx+1,1:ny,1:nz)=s2zn(1:nx+1,1:ny,1:nz,nnng)
        s3x(1:nx,1:ny+1,1:nz)=s3xn(1:nx,1:ny+1,1:nz,nnng)
        s3y(1:nx,1:ny+1,1:nz)=s3yn(1:nx,1:ny+1,1:nz,nnng)
        s3z(1:nx,1:ny+1,1:nz)=s3zn(1:nx,1:ny+1,1:nz,nnng)
        vv(1:nx,1:ny,1:nz)=vvn(1:nx,1:ny,1:nz,nnng)

        s2x(0,1:ny,1:nz)=s2x(1,1:ny,1:nz)
        s2y(0,1:ny,1:nz)=s2y(1,1:ny,1:nz)
        s2z(0,1:ny,1:nz)=s2z(1,1:ny,1:nz)
        s2x(nx+2,1:ny,1:nz)=s2x(nx+1,1:ny,1:nz)
        s2y(nx+2,1:ny,1:nz)=s2y(nx+1,1:ny,1:nz)
        s2z(nx+2,1:ny,1:nz)=s2z(nx+1,1:ny,1:nz)

        s3x(1:nx,0,1:nz)=s3x(1:nx,1,1:nz)
        s3y(1:nx,0,1:nz)=s3y(1:nx,1,1:nz)
        s3z(1:nx,0,1:nz)=s3z(1:nx,1,1:nz)
        s3x(1:nx,ny+2,1:nz)=s3x(1:nx,ny+1,1:nz)
        s3y(1:nx,ny+2,1:nz)=s3y(1:nx,ny+1,1:nz)
        s3z(1:nx,ny+2,1:nz)=s3z(1:nx,ny+1,1:nz)

        s1x(1:nx,1:ny,0)=s1x(1:nx,1:ny,1)
        s1y(1:nx,1:ny,0)=s1y(1:nx,1:ny,1)
        s1z(1:nx,1:ny,0)=s1z(1:nx,1:ny,1)
        s1x(1:nx,1:ny,nz+2)=s1x(1:nx,1:ny,nz+1)
        s1y(1:nx,1:ny,nz+2)=s1y(1:nx,1:ny,nz+1)
        s1z(1:nx,1:ny,nz+2)=s1z(1:nx,1:ny,nz+1)
    end subroutine geon
!网格参数的计算
subroutine geo1(x0,y0,z0,s01x,s01y,s01z,xx001,yy001,zz001,s02x,s02y,s02z,xx002,yy002,zz002,s03x,s03y,s03z,xx003,yy003,zz003,vv0,xx000,yy000,zz000)
    use global
    implicit none
    real(8) ::x24,y24,z24,x31,y31,z31
    real(8),intent(in)  ::x0(nx+1,ny+1,nz+1),y0(nx+1,ny+1,nz+1),z0(nx+1,ny+1,nz+1)
    real(8),intent(out) ::s01x(nx,ny,nz+1),s01y(nx,ny,nz+1),s01z(nx,ny,nz+1),xx001(nx,ny,nz+1),yy001(nx,ny,nz+1),zz001(nx,ny,nz+1),s02x(nx+1,ny,nz),s02y(nx+1,ny,nz),&
                          s02z(nx+1,ny,nz),xx002(nx+1,ny,nz),yy002(nx+1,ny,nz),zz002(nx+1,ny,nz),s03x(nx,ny+1,nz),s03y(nx,ny+1,nz),s03z(nx,ny+1,nz),xx003(nx,ny+1,nz),&
                          yy003(nx,ny+1,nz),zz003(nx,ny+1,nz),vv0(nx,ny,nz),xx000(nx,ny,nz),yy000(nx,ny,nz),zz000(nx,ny,nz)
    do k=1,nz+1
        do j=1,ny
            do i=1,nx
               x24=x0(i,j+1,k)-x0(i+1,j,k)
               y24=y0(i,j+1,k)-y0(i+1,j,k)
               z24=z0(i,j+1,k)-z0(i+1,j,k)

               x31=x0(i+1,j+1,k)-x0(i,j,k)
               y31=y0(i+1,j+1,k)-y0(i,j,k)
               z31=z0(i+1,j+1,k)-z0(i,j,k)

               s01x(i,j,k)=0.5d0*(y24*z31-z24*y31)
               s01y(i,j,k)=0.5d0*(z24*x31-x24*z31)
               s01z(i,j,k)=0.5d0*(x24*y31-y24*x31)

               xx001(i,j,k)=0.25d0*(x0(i,j,k)+x0(i,j+1,k)+x0(i+1,j,k)+x0(i+1,j+1,k))
               yy001(i,j,k)=0.25d0*(y0(i,j,k)+y0(i,j+1,k)+y0(i+1,j,k)+y0(i+1,j+1,k))
               zz001(i,j,k)=0.25d0*(z0(i,j,k)+z0(i,j+1,k)+z0(i+1,j,k)+z0(i+1,j+1,k))
            end do
        end do
    end do

    do k=1,nz
        do j=1,ny
            do i=1,nx+1
               x24=x0(i,j+1,k+1)-x0(i,j,k)
               y24=y0(i,j+1,k+1)-y0(i,j,k)
               z24=z0(i,j+1,k+1)-z0(i,j,k)

               x31=x0(i,j+1,k)-x0(i,j,k+1)
               y31=y0(i,j+1,k)-y0(i,j,k+1)
               z31=z0(i,j+1,k)-z0(i,j,k+1)

               s02x(i,j,k)=0.5d0*(y24*z31-z24*y31)
               s02y(i,j,k)=0.5d0*(z24*x31-x24*z31)
               s02z(i,j,k)=0.5d0*(x24*y31-y24*x31)

               xx002(i,j,k)=0.25d0*(x0(i,j,k)+x0(i,j+1,k)+x0(i,j,k+1)+x0(i,j+1,k+1))
               yy002(i,j,k)=0.25d0*(y0(i,j,k)+y0(i,j+1,k)+y0(i,j,k+1)+y0(i,j+1,k+1))
               zz002(i,j,k)=0.25d0*(z0(i,j,k)+z0(i,j+1,k)+z0(i,j,k+1)+z0(i,j+1,k+1))
            end do
        end do
    end do

    do k=1,nz
        do j=1,ny+1
            do i=1,nx
               x24=x0(i+1,j,k+1)-x0(i,j,k)
               y24=y0(i+1,j,k+1)-y0(i,j,k)
               z24=z0(i+1,j,k+1)-z0(i,j,k)

               x31=x0(i,j,k+1)-x0(i+1,j,k)
               y31=y0(i,j,k+1)-y0(i+1,j,k)
               z31=z0(i,j,k+1)-z0(i+1,j,k)

               s03x(i,j,k)=0.5d0*(y24*z31-z24*y31)
               s03y(i,j,k)=0.5d0*(z24*x31-x24*z31)
               s03z(i,j,k)=0.5d0*(x24*y31-y24*x31)

               xx003(i,j,k)=0.25d0*(x0(i,j,k)+x0(i+1,j,k)+x0(i,j,k+1)+x0(i+1,j,k+1))
               yy003(i,j,k)=0.25d0*(y0(i,j,k)+y0(i+1,j,k)+y0(i,j,k+1)+y0(i+1,j,k+1))
               zz003(i,j,k)=0.25d0*(z0(i,j,k)+z0(i+1,j,k)+z0(i,j,k+1)+z0(i+1,j,k+1))
            end do
        end do
    end do

    do k=1,nz
        do j=1,ny
            do i=1,nx
               x24=x0(i+1,j+1,k+1)-x0(i,j,k)
               y24=y0(i+1,j+1,k+1)-y0(i,j,k)
               z24=z0(i+1,j+1,k+1)-z0(i,j,k)

               vv0(i,j,k)=-(x24*(s01x(i,j,k)+s02x(i,j,k)+s03x(i,j,k))  &
                           +y24*(s01y(i,j,k)+s02y(i,j,k)+s03y(i,j,k))  &
                           +z24*(s01z(i,j,k)+s02z(i,j,k)+s03z(i,j,k)))/3.d0
               xx000(i,j,k)=0.125d0*(x0(i,j,k)+x0(i+1,j,k)+x0(i,j,k+1)+x0(i+1,j,k+1)+x0(i,j+1,k)+x0(i+1,j+1,k)+x0(i,j+1,k+1)+x0(i+1,j+1,k+1))
               yy000(i,j,k)=0.125d0*(y0(i,j,k)+y0(i+1,j,k)+y0(i,j,k+1)+y0(i+1,j,k+1)+y0(i,j+1,k)+y0(i+1,j+1,k)+y0(i,j+1,k+1)+y0(i+1,j+1,k+1))
               zz000(i,j,k)=0.125d0*(z0(i,j,k)+z0(i+1,j,k)+z0(i,j,k+1)+z0(i+1,j,k+1)+z0(i,j+1,k)+z0(i+1,j+1,k)+z0(i,j+1,k+1)+z0(i+1,j+1,k+1))
            end do
        end do
    end do
    end subroutine geo1

subroutine ye   !边界值
    use global
    implicit none
    integer ::kk,j1,k1
    real(8) ::vall1,vall2,vall3,vall4,val,qq2,t1,cvlt,cvl,cvu,tem,tur1,tur2,tur3,fv1

    if(nng>nng0)then !若不为粗网格
        ny=nny(nng)
        nz=nnz(nng)
        do k=1,nz
            do j=1,ny
                y1=yy2(1,j,k,nng)
                z1=zz2(1,j,k,nng)
                rr=sqrt(y1*y1+z1*z1)
                sir=z1/rr
                cor=y1/rr
                betaxn(j,k,nng,1:nt)=vxx/vee
                betayn(j,k,nng,1:nt)=(vrr*cor-vtt*sir)/vee
                betazn(j,k,nng,1:nt)=(vrr*sir+vtt*cor)/vee
                hatn(j,k,nng,1:nt)=ht
                petn(j,k,nng,1:nt)=pt
            end do
        end do
    end if
    end subroutine ye

subroutine yen(nnng) !当前网格层边界值
    use global
    implicit none
    integer,intent(in) ::nnng

        betax(1:ny,1:nz)=betaxn(1:ny,1:nz,nnng,n)
        betay(1:ny,1:nz)=betayn(1:ny,1:nz,nnng,n)
        betaz(1:ny,1:nz)=betazn(1:ny,1:nz,nnng,n)
        pet(1:ny,1:nz)=petn(1:ny,1:nz,nnng,n)
        hat(1:ny,1:nz)=hatn(1:ny,1:nz,nnng,n)
    end subroutine yen

subroutine tsd      !时间导数项计算（除湍流输运方程）
    use global
    implicit none
    integer ::nn

    ts1=0.d0
    ts2=0.d0
    ts3=0.d0
    ts4=0.d0
    ts5=0.d0
    do k=1,nz
        do j=1,ny
            do i=1,nx
                do nn=1,nt  !随着物理时间步推进，进行累加
                    ts1(i,j,k)=ts1(i,j,k)+dm(n,nn)*Q11(i,j,k,nn) !q11原始守衡量
                    ts2(i,j,k)=ts2(i,j,k)+dm(n,nn)*Q12(i,j,k,nn)
                    ts3(i,j,k)=ts3(i,j,k)+dm(n,nn)*Q13(i,j,k,nn)
                    ts4(i,j,k)=ts4(i,j,k)+dm(n,nn)*Q14(i,j,k,nn)
                    ts5(i,j,k)=ts5(i,j,k)+dm(n,nn)*Q15(i,j,k,nn)
                end do
            end do
        end do
    end do
    end subroutine tsd

subroutine tsdsa  !SA湍流输运方程（以下简称SA方程）时间导数项计算
    use global
    implicit none
    integer ::nn

    ts6=0.d0
    do k=1,nz
        do j=1,ny
            do i=1,nx
                do nn=1,nt
                    ts6(i,j,k)=ts6(i,j,k)+dm(n,nn)*Q16(i,j,k,nn)
                end do
            end do
        end do
    end do
    end subroutine tsdsa

subroutine march1 !3步R-K推进
    use global
    implicit none

    q01(1:nx,1:ny,1:nz)=q11(1:nx,1:ny,1:nz,n)
    q02(1:nx,1:ny,1:nz)=q12(1:nx,1:ny,1:nz,n)
    q03(1:nx,1:ny,1:nz)=q13(1:nx,1:ny,1:nz,n)
    q04(1:nx,1:ny,1:nz)=q14(1:nx,1:ny,1:nz,n)
    q05(1:nx,1:ny,1:nz)=q15(1:nx,1:ny,1:nz,n)
    ta=1.d0
    timl=0.6d0*cfl !时间步长
    call ppp!由守恒量计算原始量
    call bc!边界条件
    call step!当地时间步长
    call ddd!人工粘性
    call qqq!对流通量的计算
    call qqqv!黏性通量
    call pred(1)!每步R-K推进方法

    timl=0.6d0*cfl
    call ppp
    call bc
    call qqq
    call pred(0)

    timl=cfl
    call ppp
    call bc
    call qqq
    call pred(1)
    end subroutine march1

subroutine march2 !粗网格3步R-K推进
    use global
    implicit none

    q01(1:nx,1:ny,1:nz)=q11(1:nx,1:ny,1:nz,n) !R-K中存放上一次时间推进结果（推进初值）
    q02(1:nx,1:ny,1:nz)=q12(1:nx,1:ny,1:nz,n)
    q03(1:nx,1:ny,1:nz)=q13(1:nx,1:ny,1:nz,n)
    q04(1:nx,1:ny,1:nz)=q14(1:nx,1:ny,1:nz,n)
    q05(1:nx,1:ny,1:nz)=q15(1:nx,1:ny,1:nz,n)
    ta=1.d0

    timl=0.6d0*cfl
    call ppp
    call bc
    call step
    call dddc
    call qqq
    call qqqv
    call pred(1)

    timl=0.6d0*cfl
    call ppp
    call bc
    call qqq
    call pred(0)

    timl=cfl
    call ppp
    call bc
    call qqq
    call pred(1)
    end subroutine march2

subroutine marchsa!SA 3步R-K推进
    use global
    implicit none

    call tsdsa
    q06(1:nx,1:ny,1:nz)=q16(1:nx,1:ny,1:nz,n)
    ta=1.5d0

    timl=0.6d0*cfl
    call ppp
    call bc
    call step
    call dddsa
    call qqqsa
    call qqqvsa
    call predsa(1)

    timl=0.6d0*cfl
    call ppp
    call bc
    call qqqsa
    call predsa(0)

    timl=cfl
    call ppp
    call bc
    call qqqsa
    call predsa(1)

    q06(1:nx,1:ny,1:nz)=q16(1:nx,1:ny,1:nz,n)
    timl=cfl*0.125d0
    call ppp
    call bc
    call step
    call dddsa
    call qqqsa
    call qqqvsa
    call predsa(0)
    end subroutine marchsa

subroutine ppp  !求原始变量q16
    use global
    implicit none
    real(8) ::qq2,cvl

    do k=1,nz
        do j=1,ny
            do i=1,nx
                dim=Q11(i,j,k,n)
                pvx(i,j,k)=Q12(i,j,k,n)/dim !pvx/pvy/pvz原始量,q1原始守恒量
                pvy(i,j,k)=Q13(i,j,k,n)/dim
                pvz(i,j,k)=Q14(i,j,k,n)/dim
                vx=pvx(i,j,k)
                vy=pvy(i,j,k)
                vz=pvz(i,j,k)

                y1=yy0(i,j,k)
                z1=zz0(i,j,k)
                rr=sqrt(y1*y1+z1*z1)
                sir=z1/rr
                cor=y1/rr
                vth(i,j,k)=vz*cor-vy*sir !斜率，导数
                vre(i,j,k)=vz*sir+vy*cor !法向

                qq2=vx*vx+vy*vy+vz*vz !距离
                p(i,j,k)=0.4d0*(Q15(i,j,k,n)-0.5d0*dim*qq2)
                t(i,j,k)=p(i,j,k)/(dim*rg)
                cvl=cvl0*(t(i,j,k)/t0)**1.5*(t0+ts)/(t(i,j,k)+ts)
                q16(i,j,k,n)=max(q16(i,j,k,n),1.D-4*cvl)
            end do
        end do
    end do
    end subroutine ppp

subroutine bc !边界条件的设定
    use global
    implicit none
    real(8) ::qq2,sxn,syn,szn,ds,a,deltp,rrhoc
    real(8) ::u,v,w,uabs,unorm,rinv,c02,dis,cb,cosa,hb,cc02,cvl
    real(8) ::ss,sr,sv,sd,s1,ve,v2,dr,x1,y2,z2,rr22,cor2,sir2,rp1,rp2
    real(8) ::ca,qz1,qz2,qq,temp,t1,t2,cvlt,cvu,tem,tur1,tur2,tur3,fv1
    real(8),allocatable ::hr(:),hv(:),hd(:),hp(:)
    !***************进口边界***************
        do k=1,nz
            do j=1,ny
                dim=Q11(1,j,k,n) !同ppp
                vx=Q12(1,j,k,n)/dim
                vy=Q13(1,j,k,n)/dim
                vz=Q14(1,j,k,n)/dim
                qq2=vx*vx+vy*vy+vz*vz
                en=Q15(1,j,k,n)
                pp=0.4d0*(en-0.5d0*dim*qq2)
                t1=pp/(dim*rg)              !rg湍流计算所需参数
                cvlt=cvl0*(t1/t0)**1.5*(t0+ts)/(t1+ts)
                cvl=cvlt/dim
                cvu=c2*cvl
                tur1=1.D-4
                tur2=1d-6
                tur3=1d-6
                do while(abs((tur1-tur3)/tur3)>1d-6) !abs求绝对值，若满足一定的误差
                    tem=tur1/cvl
                    fv1=1.d0/(1.d0+(cv1/tem)**3)
                    tur2=cvu/fv1
                    tur3=tur1
                    tur1=tur2
                end do
                turi(j,k)=tur2

                ds=sqrt(s2x(1,j,k)*s2x(1,j,k)+s2y(1,j,k)*s2y(1,j,k)+s2z(1,j,k)*s2z(1,j,k))
                sxn=s2x(1,j,k)/ds
                syn=s2y(1,j,k)/ds
                szn=s2z(1,j,k)/ds
                u=pvx(1,j,k)
                v=pvy(1,j,k)
                w=pvz(1,j,k)
                uabs=sqrt(u*u+v*v+w*w)
                unorm = u*sxn + v*syn +w*szn
                a=sqrt(1.4d0*p(1,j,k)/q11(1,j,k,n))
                if (uabs < 1.D-20) then
                    cosa = 1.D0
                else
                    cosa = -unorm/uabs
                end if
                rinv = unorm - 5.d0*a
                c02  = a*a + 0.2d0*uabs*uabs
                dis  = (0.4d0*cosa*cosa+2.d0)*c02/(0.4d0*rinv*rinv) - 0.2d0
                if (dis < 0.D0) then
                    dis = 1.D-20
                end if
                cb    = -rinv*(0.4d0/(0.4d0*cosa*cosa+2.d0))*(1.d0+cosa*Sqrt(dis))
                cc02  = Min(cb*cb/c02, 1.d0)
                hb    = hat(j,k)*cc02

                t(0,j,k)=hb/cp
                p(0,j,k)=pet(j,k)*(cc02**3.5)
                q11(0,j,k,n)=p(0,j,k)/(t(0,j,k)*rg)
                uabs =2.d0*(hat(j,k)-hb)
                q15(0,j,k,n) =2.5d0*p(0,j,k)+0.5d0*q11(0,j,k,n)*uabs
                pvx(0,j,k)=sqrt(uabs)*betax(j,k)
                pvy(0,j,k)=sqrt(uabs)*betay(j,k)
                pvz(0,j,k)=sqrt(uabs)*betaz(j,k)
                q12(0,j,k,n)=q11(0,j,k,n)*pvx(0,j,k)
                q13(0,j,k,n)=q11(0,j,k,n)*pvy(0,j,k)
                q14(0,j,k,n)=q11(0,j,k,n)*pvz(0,j,k)
                q16(0,j,k,n)=turi(j,k)*q11(0,j,k,n)
            end do
        end do
    !***************出口边界***************
          allocate(hr(ny))
          allocate(hv(ny))
           allocate(hd(ny))
           allocate(hp(ny))
        do j=1,ny
            ss=0.d0
            sr=0.d0
            sv=0.d0
            sd=0.d0
            do k=1,nz
                s1=sqrt(s2x(nx+1,j,k)*s2x(nx+1,j,k)+s2y(nx+1,j,k)*s2y(nx+1,j,k)+s2z(nx+1,j,k)*s2z(nx+1,j,k))
                y1=yy02(nx+1,j,k)
                z1=zz02(nx+1,j,k)
                rr=sqrt(y1*y1+z1*z1)
                dim=q11(nx,j,k,n)
                ve=vth(nx,j,k)
                ss=ss+s1
                sr=sr+rr*s1
                sd=sd+dim*s1
                sv=sv+ve*s1
            end do
            hr(j)=sr/ss
            hv(j)=sv/ss
            hd(j)=sd/ss
        end do
        hp(1)=pb1
        do j=2,ny
            dim=0.5d0*(hd(j-1)+hd(j))
            v2=0.5d0*(hv(j-1)+hv(j))
            dr=hr(j)-hr(j-1)
            rr=0.5d0*(hr(j)+hr(j-1))
            hp(j)=hp(j-1)+dim*v2*v2/rr*dr
        end do
        do k=1,nz
            do j=1,ny
                peb(j,k)=hp(j)
            end do
        end do
        deallocate(hr)
        deallocate(hv)
        deallocate(hd)
        deallocate(hp)
        do k=1,nz
            do j=1,ny
                ds=sqrt(s2x(nx+1,j,k)*s2x(nx+1,j,k)+s2y(nx+1,j,k)*s2y(nx+1,j,k)+s2z(nx+1,j,k)*s2z(nx+1,j,k))
                sxn=-s2x(nx+1,j,k)/ds
                syn=-s2y(nx+1,j,k)/ds
                szn=-s2z(nx+1,j,k)/ds
                a=sqrt(1.4d0*p(nx,j,k)/q11(nx,j,k,n))
                rrhoc=1.d0/(q11(nx,j,k,n)*a)
                deltp=p(nx,j,k)-peb(j,k)

                p(nx+1,j,k)=peb(j,k)
                q11(nx+1,j,k,n) =q11(nx,j,k,n)-deltp/(a*a)
                pvx(nx+1,j,k)   =pvx(nx,j,k)+sxn*deltp*rrhoc
                pvy(nx+1,j,k)   =pvy(nx,j,k)+syn*deltp*rrhoc
                pvz(nx+1,j,k)   =pvz(nx,j,k)+szn*deltp*rrhoc
                q12(nx+1,j,k,n)=q11(nx+1,j,k,n)*pvx(nx+1,j,k)
                q13(nx+1,j,k,n)=q11(nx+1,j,k,n)*pvy(nx+1,j,k)
                q14(nx+1,j,k,n)=q11(nx+1,j,k,n)*pvz(nx+1,j,k)
                qq2=pvx(nx+1,j,k)*pvx(nx+1,j,k)+pvy(nx+1,j,k)*pvy(nx+1,j,k)+pvz(nx+1,j,k)*pvz(nx+1,j,k)
                q15(nx+1,j,k,n) =2.5d0*p(nx+1,j,k)+0.5d0*q11(nx+1,j,k,n)*qq2
                q16(nx+1,j,k,n) =q16(nx,j,k,n)
                t(nx+1,j,k)=p(nx+1,j,k)/(q11(nx+1,j,k,n)*rg)
            end do
        end do
!***************上下边界y***************
        do i=1,nx
            do k=1,nz
                !*****叶根位置****
                p(i,0,k)  = p(i,1,k)
                t(i,0,k)  = t(i,1,k)
                q11(i,0,k,n)=q11(i,1,k,n)
                pvx(i,0,k)=-pvx(i,1,k)
                pvy(i,0,k)=-2.d0*rpm*zz03(i,1,k)-pvy(i,1,k)
                pvz(i,0,k)= 2.d0*rpm*yy03(i,1,k)-pvz(i,1,k)
                q12(i,0,k,n)=q11(i,0,k,n)*pvx(i,0,k)
                q13(i,0,k,n)=q11(i,0,k,n)*pvy(i,0,k)
                q14(i,0,k,n)=q11(i,0,k,n)*pvz(i,0,k)
                qq2=pvx(i,0,k)*pvx(i,0,k)+pvy(i,0,k)*pvy(i,0,k)+pvz(i,0,k)*pvz(i,0,k)
                q15(i,0,k,n)=2.5d0*p(i,0,k)+0.5d0*q11(i,0,k,n)*qq2
                q16(i,0,k,n)=-q16(i,1,k,n)
                !*****叶顶位置****
                p(i,ny+1,k)  = p(i,ny,k)
                t(i,ny+1,k)  = t(i,ny,k)
                q11(i,ny+1,k,n)=q11(i,ny,k,n)
                pvx(i,ny+1,k)=-pvx(i,ny,k)
                pvy(i,ny+1,k)=-pvy(i,ny,k)
                pvz(i,ny+1,k)=-pvz(i,ny,k)
                q12(i,ny+1,k,n)=q11(i,ny+1,k,n)*pvx(i,ny+1,k)
                q13(i,ny+1,k,n)=q11(i,ny+1,k,n)*pvy(i,ny+1,k)
                q14(i,ny+1,k,n)=q11(i,ny+1,k,n)*pvz(i,ny+1,k)
                qq2=pvx(i,ny+1,k)*pvx(i,ny+1,k)+pvy(i,ny+1,k)*pvy(i,ny+1,k)+pvz(i,ny+1,k)*pvz(i,ny+1,k)
                q15(i,ny+1,k,n)=2.5d0*p(i,ny+1,k)+0.5d0*q11(i,ny+1,k,n)*qq2
                q16(i,ny+1,k,n)=-q16(i,ny,k,n)
            end do
        end do
!***************前后固壁边界z的计算***************
        q11(ib:it,jb:jt,0,n)  =q11(ib:it,jb:jt,1,n)
        pvx(ib:it,jb:jt,0)  =-pvx(ib:it,jb:jt,1)
        pvy(ib:it,jb:jt,0)  =-2.d0*rpm*zz01(ib:it,jb:jt,1)-pvy(ib:it,jb:jt,1)
        pvz(ib:it,jb:jt,0)  = 2.d0*rpm*yy01(ib:it,jb:jt,1)-pvz(ib:it,jb:jt,1)
        p(ib:it,jb:jt,0)  =p(ib:it,jb:jt,1)
        q16(ib:it,jb:jt,0,n)  =-q16(ib:it,jb:jt,1,n)
        q11(ib:it,jb:jt,nz+1,n)  =q11(ib:it,jb:jt,nz,n)
        pvx(ib:it,jb:jt,nz+1)  =-pvx(ib:it,jb:jt,nz)
        pvy(ib:it,jb:jt,nz+1)  =-2.d0*rpm*zz01(ib:it,jb:jt,nz+1)-pvy(ib:it,jb:jt,nz)
        pvz(ib:it,jb:jt,nz+1)  = 2.d0*rpm*yy01(ib:it,jb:jt,nz+1)-pvz(ib:it,jb:jt,nz)
        p(ib:it,jb:jt,nz+1)  =p(ib:it,jb:jt,nz)
        q16(ib:it,jb:jt,nz+1,n)  =-q16(ib:it,jb:jt,nz,n)
    !*********周期性边界上虚拟的物理量
    if(myid==0)then
        myidl=numprocs-1
        myidr=myid+1
    else if(myid==numprocs-1)then
        myidl=myid-1
        myidr=0
    else
        myidl=myid-1
        myidr=myid+1
    end if
        call MPI_SENDRECV(q11(1:ib-1,1:ny,1,n),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidl,20,q11(1:ib-1,1:ny,nz+1,n),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidr,20,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvx(1:ib-1,1:ny,1),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidl,21,pvx(1:ib-1,1:ny,nz+1),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidr,21,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(vth(1:ib-1,1:ny,1),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidl,22,vth(1:ib-1,1:ny,nz+1),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidr,22,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(vre(1:ib-1,1:ny,1),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidl,23,vre(1:ib-1,1:ny,nz+1),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidr,23,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(p(1:ib-1,1:ny,1),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidl,24,p(1:ib-1,1:ny,nz+1),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidr,24,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(q16(1:ib-1,1:ny,1,n),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidl,25,q16(1:ib-1,1:ny,nz+1,n),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidr,25,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(q11(it+1:nx,1:ny,1,n),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidl,30,q11(it+1:nx,1:ny,nz+1,n),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidr,30,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvx(it+1:nx,1:ny,1),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidl,31,pvx(it+1:nx,1:ny,nz+1),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidr,31,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(vth(it+1:nx,1:ny,1),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidl,32,vth(it+1:nx,1:ny,nz+1),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidr,32,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(vre(it+1:nx,1:ny,1),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidl,33,vre(it+1:nx,1:ny,nz+1),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidr,33,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(p(it+1:nx,1:ny,1),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidl,34,p(it+1:nx,1:ny,nz+1),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidr,34,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(q16(it+1:nx,1:ny,1,n),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidl,35,q16(it+1:nx,1:ny,nz+1,n),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidr,35,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(q11(ib:it,jt+1:ny,1,n),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidl,40,q11(ib:it,jt+1:ny,nz+1,n),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidr,40,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvx(ib:it,jt+1:ny,1),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidl,41,pvx(ib:it,jt+1:ny,nz+1),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidr,41,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(vth(ib:it,jt+1:ny,1),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidl,42,vth(ib:it,jt+1:ny,nz+1),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidr,42,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(vre(ib:it,jt+1:ny,1),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidl,43,vre(ib:it,jt+1:ny,nz+1),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidr,43,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(p(ib:it,jt+1:ny,1),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidl,44,p(ib:it,jt+1:ny,nz+1),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidr,44,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(q16(ib:it,jt+1:ny,1,n),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidl,45,q16(ib:it,jt+1:ny,nz+1,n),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidr,45,MPI_COMM_WORLD,status,ierr)

        call MPI_SENDRECV(q11(1:ib-1,1:ny,nz,n),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidr,50,q11(1:ib-1,1:ny,0,n),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidl,50,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvx(1:ib-1,1:ny,nz),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidr,51,pvx(1:ib-1,1:ny,0),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidl,51,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(vth(1:ib-1,1:ny,nz),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidr,52,vth(1:ib-1,1:ny,0),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidl,52,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(vre(1:ib-1,1:ny,nz),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidr,53,vre(1:ib-1,1:ny,0),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidl,53,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(p(1:ib-1,1:ny,nz),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidr,54,p(1:ib-1,1:ny,0),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidl,54,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(q16(1:ib-1,1:ny,nz,n),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidr,55,q16(1:ib-1,1:ny,0,n),(ib-1)*ny,MPI_DOUBLE_PRECISION,myidl,55,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(q11(it+1:nx,1:ny,nz,n),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidr,60,q11(it+1:nx,1:ny,0,n),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidl,60,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvx(it+1:nx,1:ny,nz),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidr,61,pvx(it+1:nx,1:ny,0),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidl,61,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(vth(it+1:nx,1:ny,nz),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidr,62,vth(it+1:nx,1:ny,0),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidl,62,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(vre(it+1:nx,1:ny,nz),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidr,63,vre(it+1:nx,1:ny,0),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidl,63,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(p(it+1:nx,1:ny,nz),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidr,64,p(it+1:nx,1:ny,0),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidl,64,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(q16(it+1:nx,1:ny,nz,n),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidr,65,q16(it+1:nx,1:ny,0,n),(nx-it)*ny,MPI_DOUBLE_PRECISION,myidl,65,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(q11(ib:it,jt+1:ny,nz,n),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidr,70,q11(ib:it,jt+1:ny,0,n),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidl,70,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvx(ib:it,jt+1:ny,nz),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidr,71,pvx(ib:it,jt+1:ny,0),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidl,71,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(vth(ib:it,jt+1:ny,nz),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidr,72,vth(ib:it,jt+1:ny,0),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidl,72,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(vre(ib:it,jt+1:ny,nz),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidr,73,vre(ib:it,jt+1:ny,0),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidl,73,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(p(ib:it,jt+1:ny,nz),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidr,74,p(ib:it,jt+1:ny,0),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidl,74,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(q16(ib:it,jt+1:ny,nz,n),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidr,75,q16(ib:it,jt+1:ny,0,n),(it-ib+1)*(ny-jt),MPI_DOUBLE_PRECISION,myidl,75,MPI_COMM_WORLD,status,ierr)
!*********周向z
    do j=1,ny
        do i=1,nx
            if(i>=ib .and. i<=it .and. j>=jb .and. j<=jt)then
            else
                y1=yy0(i,j,0)
                z1=zz0(i,j,0)
                rr=sqrt(y1*y1+z1*z1)
                sir=z1/rr
                cor=y1/rr
                pvy(i,j,0)=vre(i,j,0)*cor-vth(i,j,0)*sir
                pvz(i,j,0)=vre(i,j,0)*sir+vth(i,j,0)*cor

                y1=yy0(i,j,nz+1)
                z1=zz0(i,j,nz+1)
                rr=sqrt(y1*y1+z1*z1)
                sir=z1/rr
                cor=y1/rr
                pvy(i,j,nz+1)=vre(i,j,nz+1)*cor-vth(i,j,nz+1)*sir
                pvz(i,j,nz+1)=vre(i,j,nz+1)*sir+vth(i,j,nz+1)*cor
            end if
            q12(i,j,0,n)= q11(i,j,0,n)*pvx(i,j,0)
            q13(i,j,0,n)= q11(i,j,0,n)*pvy(i,j,0)
            q14(i,j,0,n)= q11(i,j,0,n)*pvz(i,j,0)
            qq2=pvx(i,j,0)*pvx(i,j,0)+pvy(i,j,0)*pvy(i,j,0)+pvz(i,j,0)*pvz(i,j,0)
            q15(i,j,0,n)=2.5d0*p(i,j,0)+0.5d0*q11(i,j,0,n)*qq2
            t(i,j,0)=p(i,j,0)/(Q11(i,j,0,n)*rg)

            q12(i,j,nz+1,n)= q11(i,j,nz+1,n)*pvx(i,j,nz+1)
            q13(i,j,nz+1,n)= q11(i,j,nz+1,n)*pvy(i,j,nz+1)
            q14(i,j,nz+1,n)= q11(i,j,nz+1,n)*pvz(i,j,nz+1)
            qq2=pvx(i,j,nz+1)*pvx(i,j,nz+1)+pvy(i,j,nz+1)*pvy(i,j,nz+1)+pvz(i,j,nz+1)*pvz(i,j,nz+1)
            q15(i,j,nz+1,n)=2.5d0*p(i,j,nz+1)+0.5d0*q11(i,j,nz+1,n)*qq2
            t(i,j,nz+1)=p(i,j,nz+1)/(Q11(i,j,nz+1,n)*rg)
        end do
    end do
    end subroutine bc

subroutine step !当地时间步长 ，同时求解谱半径
    use global
    implicit none
    real(8) :: tc,td,aa,cv,kc
    real(8),allocatable ::sri1(:,:,:),srj1(:,:,:),srk1(:,:,:)

    do k=1,nz
        do j=1,ny
            do i=1,nx
                !******x方向
                vx=0.5d0*(pvx(i,j,k)+pvx(i+1,j,k))
                vy=0.5d0*(pvy(i,j,k)+pvy(i+1,j,k))
                vz=0.5d0*(pvz(i,j,k)+pvz(i+1,j,k))
                wx=vx      !只用到了vx？？？？？？
                wy=vy+rpm*zz02(i+1,j,k)
                wz=vz-rpm*yy02(i+1,j,k)
                tc=abs(wx*s2x(i+1,j,k)+wy*s2y(i+1,j,k)+wz*s2z(i+1,j,k))
                !******y方向
                vx=0.5d0*(pvx(i,j,k)+pvx(i,j+1,k))
                vy=0.5d0*(pvy(i,j,k)+pvy(i,j+1,k))
                vz=0.5d0*(pvz(i,j,k)+pvz(i,j+1,k))
                wx=vx
                wy=vy+rpm*zz03(i,j+1,k)
                wz=vz-rpm*yy03(i,j+1,k)
                tc=abs(wx*s3x(i,j+1,k)+wy*s3y(i,j+1,k)+wz*s3z(i,j+1,k))+tc
                !******z方向
                vx=0.5d0*(pvx(i,j,k)+pvx(i,j,k+1))
                vy=0.5d0*(pvy(i,j,k)+pvy(i,j,k+1))
                vz=0.5d0*(pvz(i,j,k)+pvz(i,j,k+1))
                wx=vx
                wy=vy+rpm*zz01(i,j,k+1)
                wz=vz-rpm*yy01(i,j,k+1)
                tc=abs(wx*s1x(i,j,k+1)+wy*s1y(i,j,k+1)+wz*s1z(i,j,k+1))+tc

                aa=sqrt(1.4d0*p(i,j,k)/q11(i,j,k,n))
                tc=tc+aa*(sqrt(s1x(i,j,k+1)**2+s1y(i,j,k+1)**2+s1z(i,j,k+1)**2)+sqrt(s2x(i+1,j,k)**2+&
                   s2y(i+1,j,k)**2+s2z(i+1,j,k)**2)+sqrt(s3x(i,j+1,k)**2+s3y(i,j+1,k)**2+s3z(i,j+1,k)**2))

                td=s1x(i,j,k+1)**2+s1y(i,j,k+1)**2+s1z(i,j,k+1)**2+s2x(i+1,j,k)**2+s2y(i+1,j,k)**2+&
                   s2z(i+1,j,k)**2+s3x(i,j+1,k)**2+s3y(i,j+1,k)**2+s3z(i,j+1,k)**2+2.d0*(abs(s1x(i,j,k+1)&
                   *s2x(i+1,j,k)+s1y(i,j,k+1)*s2y(i+1,j,k)+s1z(i,j,k+1)*s2z(i+1,j,k))+abs(s1x(i,j,k+1)&
                   *s3x(i,j+1,k)+s1y(i,j,k+1)*s3y(i,j+1,k)+s1z(i,j,k+1)*s3z(i,j+1,k))+abs(s3x(i,j+1,k)&
                   *s2x(i+1,j,k)+s3y(i,j+1,k)*s2y(i+1,j,k)+s3z(i,j+1,k)*s2z(i+1,j,k)))
                call viscosity(t(i,j,k),q16(i,j,k,n),cv,kc)
                td=8.d0*cv*td/q11(i,j,k,n)/vv(i,j,k)
                time(i,j,k)=tc+td
                sri(i,j,k)=tc+td
                srj(i,j,k)=tc+td
                srk(i,j,k)=tc+td
            end do
        end do
    end do
    allocate(sri1(0:nx+1,0:ny+1,0:nz+1))
    allocate(srj1(0:nx+1,0:ny+1,0:nz+1))
    allocate(srk1(0:nx+1,0:ny+1,0:nz+1))
    do k=1,nz
        do j=1,ny
            do i=1,nx
                sri1(i,j,k)=sri(i,j,k)*(1.d0+(srj(i,j,k)/sri(i,j,k))**0.5+(srk(i,j,k)/sri(i,j,k))**0.5)
                srj1(i,j,k)=srj(i,j,k)*(1.d0+(sri(i,j,k)/srj(i,j,k))**0.5+(srk(i,j,k)/srj(i,j,k))**0.5)
                srk1(i,j,k)=srk(i,j,k)*(1.d0+(sri(i,j,k)/srk(i,j,k))**0.5+(srj(i,j,k)/srk(i,j,k))**0.5)
            end do
        end do
    end do
    do k=1,nz
        do j=1,ny
            sri1(0,j,k) = sri1(1,j,k)
            srj1(0,j,k) = srj1(1,j,k)
            srk1(0,j,k) = srk1(1,j,k)

            sri1(nx+1,j,k) = sri1(nx,j,k)
            srj1(nx+1,j,k) = srj1(nx,j,k)
            srk1(nx+1,j,k) = srk1(nx,j,k)
        end do
    end do
    do k=1,nz
        do i=1,nx
            sri1(i,0,k) = sri1(i,1,k)
            srj1(i,0,k) = srj1(i,1,k)
            srk1(i,0,k) = srk1(i,1,k)

            sri1(i,ny+1,k) = sri1(i,ny,k)
            srj1(i,ny+1,k) = srj1(i,ny,k)
            srk1(i,ny+1,k) = srk1(i,ny,k)
        end do
    end do
    do j=1,ny
        do i=1,nx
            sri1(i,j,0) = sri1(i,j,1)
            srj1(i,j,0) = srj1(i,j,1)
            srk1(i,j,0) = srk1(i,j,1)

            sri1(i,j,nz+1) = sri1(i,j,nz)
            srj1(i,j,nz+1) = srj1(i,j,nz)
            srk1(i,j,nz+1) = srk1(i,j,nz)
        end do
    end do
    do k=0,nz+1
        do j=0,ny+1
            do i=0,nx+1
                sri(i,j,k) = sri1(i,j,k)
                srj(i,j,k) = srj1(i,j,k)
                srk(i,j,k) = srk1(i,j,k)
            end do
        end do
    end do
    deallocate(sri1)
    deallocate(srj1)
    deallocate(srk1)
    end subroutine step

subroutine ddd !人工粘性
    use global
    implicit none
    real(8) ::em2,em4
    real(8),allocatable ::dp(:)
    real(8) ::flu1,flu2,flu3,flu4,flu5,ram
    !av1~6人工粘性
    av1=0.d0
    av2=0.d0
    av3=0.d0
    av4=0.d0
    av5=0.d0
    q15(0:nx+1,0:ny+1,0:nz+1,n)=q15(0:nx+1,0:ny+1,0:nz+1,n)+p(0:nx+1,0:ny+1,0:nz+1)
    ! i-direction -----------------------------------------------------------------
    allocate(dp(0:nx+1))
    dp(0)=0.d0
    dp(nx+1)=0.d0
    do k=1,nz
      do j=1,ny
         do i=1,nx
             dp(i)=abs((p(i+1,j,k)-2.d0*p(i,j,k)+p(i-1,j,k))/(p(i+1,j,k)+2.d0*p(i,j,k)+p(i-1,j,k)))
         end do
        do i=0,nx
           ram   = 0.5d0*(sri(i,j,k)+sri(i+1,j,k))
           em2   =a2*ram*max(dp(i),dp(i+1))
           em4   =a4*ram
           em4   =max(em4-em2,0.d0)

           flu1  =em2*(Q11(i+1,j,k,n)-Q11(i,j,k,n))
           flu2  =em2*(Q12(i+1,j,k,n)-Q12(i,j,k,n))
           flu3  =em2*(Q13(i+1,j,k,n)-Q13(i,j,k,n))
           flu4  =em2*(Q14(i+1,j,k,n)-Q14(i,j,k,n))
           flu5  =em2*(Q15(i+1,j,k,n)-Q15(i,j,k,n))

           if(i>=1 .and. i<=nx-1) then
               flu1=flu1+em4*(q11(i-1,j,k,n)-3.d0*q11(i,j,k,n)+3.d0*q11(i+1,j,k,n)-q11(i+2,j,k,n))
               flu2=flu2+em4*(q12(i-1,j,k,n)-3.d0*q12(i,j,k,n)+3.d0*q12(i+1,j,k,n)-q12(i+2,j,k,n))
               flu3=flu3+em4*(q13(i-1,j,k,n)-3.d0*q13(i,j,k,n)+3.d0*q13(i+1,j,k,n)-q13(i+2,j,k,n))
               flu4=flu4+em4*(q14(i-1,j,k,n)-3.d0*q14(i,j,k,n)+3.d0*q14(i+1,j,k,n)-q14(i+2,j,k,n))
               flu5=flu5+em4*(q15(i-1,j,k,n)-3.d0*q15(i,j,k,n)+3.d0*q15(i+1,j,k,n)-q15(i+2,j,k,n))
           end if
           av1(i  ,j,k) = av1(i  ,j,k) +flu1
           av2(i  ,j,k) = av2(i  ,j,k) +flu2
           av3(i  ,j,k) = av3(i  ,j,k) +flu3
           av4(i  ,j,k) = av4(i  ,j,k) +flu4
           av5(i  ,j,k) = av5(i  ,j,k) +flu5

           av1(i+1,j,k) = av1(i+1,j,k) -flu1
           av2(i+1,j,k) = av2(i+1,j,k) -flu2
           av3(i+1,j,k) = av3(i+1,j,k) -flu3
           av4(i+1,j,k) = av4(i+1,j,k) -flu4
           av5(i+1,j,k) = av5(i+1,j,k) -flu5
        end do
      end do
    end do
    deallocate(dp)
    ! j-direction -----------------------------------------------------------------
    allocate(dp(0:ny+1))
    dp(0)=0.d0
    dp(ny+1)=0.d0
    do k=1,nz
      do i=1,nx
         do j=1,ny
             dp(j)=abs((p(i,j+1,k)-2.d0*p(i,j,k)+p(i,j-1,k))/(p(i,j+1,k)+2.d0*p(i,j,k)+p(i,j-1,k)))
         end do
        do j=0,ny
           ram   = 0.5d0*(srj(i,j,k)+srj(i,j+1,k))
           em2   =a2*ram*max(dp(j),dp(j+1))
           em4   =a4*ram
           em4   =max(em4-em2,0.d0)

           flu1  =em2*(Q11(i,j+1,k,n)-Q11(i,j,k,n))
           flu2  =em2*(Q12(i,j+1,k,n)-Q12(i,j,k,n))
           flu3  =em2*(Q13(i,j+1,k,n)-Q13(i,j,k,n))
           flu4  =em2*(Q14(i,j+1,k,n)-Q14(i,j,k,n))
           flu5  =em2*(Q15(i,j+1,k,n)-Q15(i,j,k,n))

           if(j>=1 .and. j<=ny-1) then
               flu1=flu1+em4*(q11(i,j-1,k,n)-3.d0*q11(i,j,k,n)+3.d0*q11(i,j+1,k,n)-q11(i,j+2,k,n))
               flu2=flu2+em4*(q12(i,j-1,k,n)-3.d0*q12(i,j,k,n)+3.d0*q12(i,j+1,k,n)-q12(i,j+2,k,n))
               flu3=flu3+em4*(q13(i,j-1,k,n)-3.d0*q13(i,j,k,n)+3.d0*q13(i,j+1,k,n)-q13(i,j+2,k,n))
               flu4=flu4+em4*(q14(i,j-1,k,n)-3.d0*q14(i,j,k,n)+3.d0*q14(i,j+1,k,n)-q14(i,j+2,k,n))
               flu5=flu5+em4*(q15(i,j-1,k,n)-3.d0*q15(i,j,k,n)+3.d0*q15(i,j+1,k,n)-q15(i,j+2,k,n))
           end if
           av1(i  ,j,k) = av1(i  ,j,k) +flu1
           av2(i  ,j,k) = av2(i  ,j,k) +flu2
           av3(i  ,j,k) = av3(i  ,j,k) +flu3
           av4(i  ,j,k) = av4(i  ,j,k) +flu4
           av5(i  ,j,k) = av5(i  ,j,k) +flu5

           av1(i,j+1,k) = av1(i,j+1,k) -flu1
           av2(i,j+1,k) = av2(i,j+1,k) -flu2
           av3(i,j+1,k) = av3(i,j+1,k) -flu3
           av4(i,j+1,k) = av4(i,j+1,k) -flu4
           av5(i,j+1,k) = av5(i,j+1,k) -flu5
        end do
      end do
    end do
    deallocate(dp)
 ! k-direction -----------------------------------------------------------------
    allocate(dp(0:nz+1))
    dp(0)=0.d0
    dp(nz+1)=0.d0
    do j=1,ny
       do i=1,nx
              do k=1,nz
                  dp(k)=abs((p(i,j,k+1)-2.d0*p(i,j,k)+p(i,j,k-1))/(p(i,j,k+1)+2.d0*p(i,j,k)+p(i,j,k-1)))
              end do
              do k=0,nz
                  ram   = 0.5d0*(srk(i,j,k)+srk(i,j,k+1))
                  em2   =a2*ram*max(dp(k),dp(k+1))
                  em4   =a4*ram
                  em4   =max(em4-em2,0.d0)

                  flu1  =em2*(Q11(i,j,k+1,n)-Q11(i,j,k,n))
                  flu2  =em2*(Q12(i,j,k+1,n)-Q12(i,j,k,n))
                  flu3  =em2*(Q13(i,j,k+1,n)-Q13(i,j,k,n))
                  flu4  =em2*(Q14(i,j,k+1,n)-Q14(i,j,k,n))
                  flu5  =em2*(Q15(i,j,k+1,n)-Q15(i,j,k,n))

                  if(k>=1 .and. k<=nz-1) then
                      flu1=flu1+em4*(q11(i,j,k-1,n)-3.d0*q11(i,j,k,n)+3.d0*q11(i,j,k+1,n)-q11(i,j,k+2,n))
                      flu2=flu2+em4*(q12(i,j,k-1,n)-3.d0*q12(i,j,k,n)+3.d0*q12(i,j,k+1,n)-q12(i,j,k+2,n))
                      flu3=flu3+em4*(q13(i,j,k-1,n)-3.d0*q13(i,j,k,n)+3.d0*q13(i,j,k+1,n)-q13(i,j,k+2,n))
                      flu4=flu4+em4*(q14(i,j,k-1,n)-3.d0*q14(i,j,k,n)+3.d0*q14(i,j,k+1,n)-q14(i,j,k+2,n))
                      flu5=flu5+em4*(q15(i,j,k-1,n)-3.d0*q15(i,j,k,n)+3.d0*q15(i,j,k+1,n)-q15(i,j,k+2,n))
                  end if
                  av1(i  ,j,k) = av1(i  ,j,k) +flu1
                  av2(i  ,j,k) = av2(i  ,j,k) +flu2
                  av3(i  ,j,k) = av3(i  ,j,k) +flu3
                  av4(i  ,j,k) = av4(i  ,j,k) +flu4
                  av5(i  ,j,k) = av5(i  ,j,k) +flu5

                  av1(i,j,k+1) = av1(i,j,k+1) -flu1
                  av2(i,j,k+1) = av2(i,j,k+1) -flu2
                  av3(i,j,k+1) = av3(i,j,k+1) -flu3
                  av4(i,j,k+1) = av4(i,j,k+1) -flu4
                  av5(i,j,k+1) = av5(i,j,k+1) -flu5
               end do
       end do
    end do
    q15(0:nx+1,0:ny+1,0:nz+1,n)=q15(0:nx+1,0:ny+1,0:nz+1,n)-p(0:nx+1,0:ny+1,0:nz+1)
    end subroutine ddd

    subroutine dddsa  !!人工粘性SA
    use global
    implicit none
    real(8) ::em2,em4
    real(8),allocatable ::dp(:)
    real(8) ::flu6,ram

    av6=0.d0
    ! i-direction -----------------------------------------------------------------
    allocate(dp(0:nx+1))
    dp(0)=0.d0
    dp(nx+1)=0.d0
    do k=1,nz
      do j=1,ny
        do i=1,nx
            dp(i)=abs((p(i+1,j,k)-2.d0*p(i,j,k)+p(i-1,j,k))/(p(i+1,j,k)+2.d0*p(i,j,k)+p(i-1,j,k)))
        end do
        do i=0,nx
           ram   = 0.5d0*(sri(i,j,k)+sri(i+1,j,k))
           em2   =a2*ram*max(dp(i),dp(i+1))
           em4   =a4*ram
           em4   =max(em4-em2,0.d0)
           flu6  =em2*(Q16(i+1,j,k,n)-Q16(i,j,k,n))
           if(i>=1 .and. i<=nx-1) then
               flu6=flu6+em4*(q16(i-1,j,k,n)-3.*q16(i,j,k,n)+3.*q16(i+1,j,k,n)-q16(i+2,j,k,n))
           end if
           av6(i  ,j,k) = av6(i  ,j,k) +flu6
           av6(i+1,j,k) = av6(i+1,j,k) -flu6
        end do
      end do
    end do
    deallocate(dp)
    ! j-direction -----------------------------------------------------------------
    allocate(dp(0:ny+1))
    dp(0)=0.d0
    dp(ny+1)=0.d0
    do k=1,nz
      do i=1,nx
         do j=1,ny
             dp(j)=abs((p(i,j+1,k)-2.d0*p(i,j,k)+p(i,j-1,k))/(p(i,j+1,k)+2.d0*p(i,j,k)+p(i,j-1,k)))
         end do
         do j=0,ny
           ram   = 0.5d0*(srj(i,j,k)+srj(i,j+1,k))
           em2   =a2*ram*max(dp(j),dp(j+1))
           em4   =a4*ram
           em4   =max(em4-em2,0.d0)
           flu6  =em2*(Q16(i,j+1,k,n)-Q16(i,j,k,n))
           if(j>=1 .and. j<=ny-1) then
               flu6=flu6+em4*(q16(i,j-1,k,n)-3.*q16(i,j,k,n)+3.*q16(i,j+1,k,n)-q16(i,j+2,k,n))
           end if
           av6(i  ,j,k) = av6(i  ,j,k) +flu6
           av6(i,j+1,k) = av6(i,j+1,k) -flu6
         end do
      end do
    end do
    deallocate(dp)
    ! k-direction -----------------------------------------------------------------
    allocate(dp(0:nz+1))
    dp(0)=0.d0
    dp(nz+1)=0.d0
    do j=1,ny
       do i=1,nx
              do k=1,nz
                  dp(k)=abs((p(i,j,k+1)-2.d0*p(i,j,k)+p(i,j,k-1))/(p(i,j,k+1)+2.d0*p(i,j,k)+p(i,j,k-1)))
              end do
              do k=0,nz
                  ram   = 0.5d0*(srk(i,j,k)+srk(i,j,k+1))
                  em2   =a2*ram*max(dp(k),dp(k+1))
                  em4   =a4*ram
                  em4   =max(em4-em2,0.d0)
                  flu6  =em2*(Q16(i,j,k+1,n)-Q16(i,j,k,n))
                  if(k>=1 .and. k<=nz-1) then
                      flu6=flu6+em4*(q16(i,j,k-1,n)-3.*q16(i,j,k,n)+3.*q16(i,j,k+1,n)-q16(i,j,k+2,n))
                  end if
                  av6(i  ,j,k) = av6(i  ,j,k) +flu6
                  av6(i,j,k+1) = av6(i,j,k+1) -flu6
               end do
       end do
    end do
    deallocate(dp)
    end subroutine dddsa

   subroutine dddc !粗网格人工粘性
    use global
    implicit none
    real(8) ::em2
    real(8) ::flu1,flu2,flu3,flu4,flu5,ram

    av1=0.d0
    av2=0.d0
    av3=0.d0
    av4=0.d0
    av5=0.d0
    q15(0:nx+1,0:ny+1,0:nz+1,n)=q15(0:nx+1,0:ny+1,0:nz+1,n)+p(0:nx+1,0:ny+1,0:nz+1)
    ! i-direction -----------------------------------------------------------------
    em2=a2/1024.d0*sri(1,1,1)
    do k=1,nz
        do j=1,ny
            do i=0,nx
                flu1  =em2*(Q11(i+1,j,k,n)-Q11(i,j,k,n))
                flu2  =em2*(Q12(i+1,j,k,n)-Q12(i,j,k,n))
                flu3  =em2*(Q13(i+1,j,k,n)-Q13(i,j,k,n))
                flu4  =em2*(Q14(i+1,j,k,n)-Q14(i,j,k,n))
                flu5  =em2*(Q15(i+1,j,k,n)-Q15(i,j,k,n))
                av1(i  ,j,k) = av1(i  ,j,k) +flu1
                av2(i  ,j,k) = av2(i  ,j,k) +flu2
                av3(i  ,j,k) = av3(i  ,j,k) +flu3
                av4(i  ,j,k) = av4(i  ,j,k) +flu4
                av5(i  ,j,k) = av5(i  ,j,k) +flu5

                av1(i+1,j,k) = av1(i+1,j,k) -flu1
                av2(i+1,j,k) = av2(i+1,j,k) -flu2
                av3(i+1,j,k) = av3(i+1,j,k) -flu3
                av4(i+1,j,k) = av4(i+1,j,k) -flu4
                av5(i+1,j,k) = av5(i+1,j,k) -flu5
            end do
        end do
    end do
    ! j-direction -----------------------------------------------------------------
    em2=a2/1024.d0*srj(1,1,1)
    do k=1,nz
        do i=1,nx
            do j=0,ny
                flu1  =em2*(Q11(i,j+1,k,n)-Q11(i,j,k,n))
                flu2  =em2*(Q12(i,j+1,k,n)-Q12(i,j,k,n))
                flu3  =em2*(Q13(i,j+1,k,n)-Q13(i,j,k,n))
                flu4  =em2*(Q14(i,j+1,k,n)-Q14(i,j,k,n))
                flu5  =em2*(Q15(i,j+1,k,n)-Q15(i,j,k,n))
                av1(i  ,j,k) = av1(i  ,j,k) +flu1
                av2(i  ,j,k) = av2(i  ,j,k) +flu2
                av3(i  ,j,k) = av3(i  ,j,k) +flu3
                av4(i  ,j,k) = av4(i  ,j,k) +flu4
                av5(i  ,j,k) = av5(i  ,j,k) +flu5

                av1(i,j+1,k) = av1(i,j+1,k) -flu1
                av2(i,j+1,k) = av2(i,j+1,k) -flu2
                av3(i,j+1,k) = av3(i,j+1,k) -flu3
                av4(i,j+1,k) = av4(i,j+1,k) -flu4
                av5(i,j+1,k) = av5(i,j+1,k) -flu5
            end do
        end do
    end do
    ! k-direction -----------------------------------------------------------------
    em2=a2/1024.d0*srk(1,1,1)
    do j=1,ny
       do i=1,nx
              do k=0,nz
                  flu1  =em2*(Q11(i,j,k+1,n)-Q11(i,j,k,n))
                  flu2  =em2*(Q12(i,j,k+1,n)-Q12(i,j,k,n))
                  flu3  =em2*(Q13(i,j,k+1,n)-Q13(i,j,k,n))
                  flu4  =em2*(Q14(i,j,k+1,n)-Q14(i,j,k,n))
                  flu5  =em2*(Q15(i,j,k+1,n)-Q15(i,j,k,n))
                  av1(i  ,j,k) = av1(i  ,j,k) +flu1
                  av2(i  ,j,k) = av2(i  ,j,k) +flu2
                  av3(i  ,j,k) = av3(i  ,j,k) +flu3
                  av4(i  ,j,k) = av4(i  ,j,k) +flu4
                  av5(i  ,j,k) = av5(i  ,j,k) +flu5

                  av1(i,j,k+1) = av1(i,j,k+1) -flu1
                  av2(i,j,k+1) = av2(i,j,k+1) -flu2
                  av3(i,j,k+1) = av3(i,j,k+1) -flu3
                  av4(i,j,k+1) = av4(i,j,k+1) -flu4
                  av5(i,j,k+1) = av5(i,j,k+1) -flu5
               end do
       end do
    end do
    q15(0:nx+1,0:ny+1,0:nz+1,n)=q15(0:nx+1,0:ny+1,0:nz+1,n)-p(0:nx+1,0:ny+1,0:nz+1)
    end subroutine dddc

subroutine qqq!对流通量的计算
    use global
    implicit none
    real(8) ::flu1,flu2,flu3,flu4,flu5,vf,rf,qq2

    qc1=0.d0
    qc2=0.d0
    qc3=0.d0
    qc4=0.d0
    qc5=0.d0
    !******************x方向********************
    do k=1,nz
        do j=1,ny
            do i=1,nx+1
               vx=0.5d0*(pvx(i,j,k)+pvx(i-1,j,k))
               vy=0.5d0*(pvy(i,j,k)+pvy(i-1,j,k))
               vz=0.5d0*(pvz(i,j,k)+pvz(i-1,j,k))
               vf=vx*s2x(i,j,k)+vy*s2y(i,j,k)+vz*s2z(i,j,k)
               rf=rpm*zz02(i,j,k)*s2y(i,j,k)-rpm*yy02(i,j,k)*s2z(i,j,k)
               vf=vf+rf
               dim=0.5d0*(Q11(i,j,k,n)+Q11(i-1,j,k,n))
               pp=0.5d0*(p(i,j,k)+p(i-1,j,k))
               en=0.5d0*(Q15(i,j,k,n)+Q15(i-1,j,k,n))

               flu1=dim*vf
               flu2=flu1*vx+pp*s2x(i,j,k)
               flu3=flu1*vy+pp*s2y(i,j,k)
               flu4=flu1*vz+pp*s2z(i,j,k)
               flu5=(en+pp)*vf-pp*rf

               qc1(i,j,k)=qc1(i,j,k)+flu1
               qc2(i,j,k)=qc2(i,j,k)+flu2
               qc3(i,j,k)=qc3(i,j,k)+flu3
               qc4(i,j,k)=qc4(i,j,k)+flu4
               qc5(i,j,k)=qc5(i,j,k)+flu5

               qc1(i-1,j,k)=qc1(i-1,j,k)-flu1
               qc2(i-1,j,k)=qc2(i-1,j,k)-flu2
               qc3(i-1,j,k)=qc3(i-1,j,k)-flu3
               qc4(i-1,j,k)=qc4(i-1,j,k)-flu4
               qc5(i-1,j,k)=qc5(i-1,j,k)-flu5
            end do
        end do
    end do
    !******************y方向********************
    do k=1,nz
        do j=1,ny+1
            do i=1,nx
               vx=0.5d0*(pvx(i,j,k)+pvx(i,j-1,k))
               vy=0.5d0*(pvy(i,j,k)+pvy(i,j-1,k))
               vz=0.5d0*(pvz(i,j,k)+pvz(i,j-1,k))
               vf=vx*s3x(i,j,k)+vy*s3y(i,j,k)+vz*s3z(i,j,k)
               rf=rpm*zz03(i,j,k)*s3y(i,j,k)-rpm*yy03(i,j,k)*s3z(i,j,k)
               vf=vf+rf
               dim=0.5d0*(Q11(i,j,k,n)+Q11(i,j-1,k,n))
               pp=0.5d0*(p(i,j,k)+p(i,j-1,k))
               en=0.5d0*(Q15(i,j,k,n)+Q15(i,j-1,k,n))

               flu1=dim*vf
               flu2=flu1*vx+pp*s3x(i,j,k)
               flu3=flu1*vy+pp*s3y(i,j,k)
               flu4=flu1*vz+pp*s3z(i,j,k)
               flu5=(en+pp)*vf-pp*rf

               qc1(i,j,k)=qc1(i,j,k)+flu1
               qc2(i,j,k)=qc2(i,j,k)+flu2
               qc3(i,j,k)=qc3(i,j,k)+flu3
               qc4(i,j,k)=qc4(i,j,k)+flu4
               qc5(i,j,k)=qc5(i,j,k)+flu5

               qc1(i,j-1,k)=qc1(i,j-1,k)-flu1
               qc2(i,j-1,k)=qc2(i,j-1,k)-flu2
               qc3(i,j-1,k)=qc3(i,j-1,k)-flu3
               qc4(i,j-1,k)=qc4(i,j-1,k)-flu4
               qc5(i,j-1,k)=qc5(i,j-1,k)-flu5
            end do
        end do
    end do
    !******************z方向********************
    do k=1,nz+1
        do j=1,ny
            do i=1,nx
                vx=0.5d0*(pvx(i,j,k)+pvx(i,j,k-1))
                vy=0.5d0*(pvy(i,j,k)+pvy(i,j,k-1))
                vz=0.5d0*(pvz(i,j,k)+pvz(i,j,k-1))
                vf=vx*s1x(i,j,k)+vy*s1y(i,j,k)+vz*s1z(i,j,k)
                rf=rpm*zz01(i,j,k)*s1y(i,j,k)-rpm*yy01(i,j,k)*s1z(i,j,k)
                vf=vf+rf
                dim=0.5d0*(Q11(i,j,k,n)+Q11(i,j,k-1,n))
                pp=0.5d0*(p(i,j,k)+p(i,j,k-1))
                en=0.5d0*(Q15(i,j,k,n)+Q15(i,j,k-1,n))

                flu1=dim*vf
                flu2=flu1*vx+pp*s1x(i,j,k)
                flu3=flu1*vy+pp*s1y(i,j,k)
                flu4=flu1*vz+pp*s1z(i,j,k)
                flu5=(en+pp)*vf-pp*rf

                qc1(i,j,k)=qc1(i,j,k)+flu1
                qc2(i,j,k)=qc2(i,j,k)+flu2
                qc3(i,j,k)=qc3(i,j,k)+flu3
                qc4(i,j,k)=qc4(i,j,k)+flu4
                qc5(i,j,k)=qc5(i,j,k)+flu5

                qc1(i,j,k-1)=qc1(i,j,k-1)-flu1
                qc2(i,j,k-1)=qc2(i,j,k-1)-flu2
                qc3(i,j,k-1)=qc3(i,j,k-1)-flu3
                qc4(i,j,k-1)=qc4(i,j,k-1)-flu4
                qc5(i,j,k-1)=qc5(i,j,k-1)-flu5
            end do
        end do
    end do
    end subroutine qqq

subroutine qqqsa!对流通量，SA
    use global
    implicit none
    real(8) ::flu6,tur,vf,rf,qq2

    qc6=0.d0
    !******************x方向********************
    do k=1,nz
        do j=1,ny
            do i=1,nx+1
               vx=0.5d0*(pvx(i,j,k)+pvx(i-1,j,k))
               vy=0.5d0*(pvy(i,j,k)+pvy(i-1,j,k))
               vz=0.5d0*(pvz(i,j,k)+pvz(i-1,j,k))
               vf=vx*s2x(i,j,k)+vy*s2y(i,j,k)+vz*s2z(i,j,k)
               rf=rpm*zz02(i,j,k)*s2y(i,j,k)-rpm*yy02(i,j,k)*s2z(i,j,k)
               vf=vf+rf
               tur=0.5d0*(Q16(i,j,k,n)+Q16(i-1,j,k,n))
               flu6=tur*vf
               qc6(i,j,k)=qc6(i,j,k)+flu6
               qc6(i-1,j,k)=qc6(i-1,j,k)-flu6
            end do
        end do
    end do
   !******************y方向********************
    do k=1,nz
        do j=1,ny+1
            do i=1,nx
               vx=0.5d0*(pvx(i,j,k)+pvx(i,j-1,k))
               vy=0.5d0*(pvy(i,j,k)+pvy(i,j-1,k))
               vz=0.5d0*(pvz(i,j,k)+pvz(i,j-1,k))
               vf=vx*s3x(i,j,k)+vy*s3y(i,j,k)+vz*s3z(i,j,k)
               rf=rpm*zz03(i,j,k)*s3y(i,j,k)-rpm*yy03(i,j,k)*s3z(i,j,k)
               vf=vf+rf
               tur=0.5d0*(Q16(i,j,k,n)+Q16(i,j-1,k,n))
               flu6=tur*vf
               qc6(i,j,k)=qc6(i,j,k)+flu6
               qc6(i,j-1,k)=qc6(i,j-1,k)-flu6
            end do
        end do
    end do
 !******************z方向********************
    do k=1,nz+1
        do j=1,ny
            do i=1,nx
                vx=0.5d0*(pvx(i,j,k)+pvx(i,j,k-1))
                vy=0.5d0*(pvy(i,j,k)+pvy(i,j,k-1))
                vz=0.5d0*(pvz(i,j,k)+pvz(i,j,k-1))
                vf=vx*s1x(i,j,k)+vy*s1y(i,j,k)+vz*s1z(i,j,k)
                rf=rpm*zz01(i,j,k)*s1y(i,j,k)-rpm*yy01(i,j,k)*s1z(i,j,k)
                vf=vf+rf
                tur=0.5d0*(Q16(i,j,k,n)+Q16(i,j,k-1,n))
                flu6=tur*vf
                qc6(i,j,k)=qc6(i,j,k)+flu6
                qc6(i,j,k-1)=qc6(i,j,k-1)-flu6
            end do
        end do
    end do
    end subroutine qqqsa

subroutine  qqqv !粘性通量计算
    use global
    implicit none
    real(8) ::flu2,flu3,flu4,flu5,tauxx,tauyy,tauzz,tauxy,tauxz,tauyz,phix,phiy,phiz
    real(8) ::two3,uav,vav,wav,q16av,tav,mav,kav

    two3=2.D0/3.D0
    call gradsface
    qv2=0.d0
    qv3=0.d0
    qv4=0.d0
    qv5=0.d0
    ! i-direction -----------------------------------------------------------------
    do k=1,nz
        do j=1,ny
            do i=1,nx+1
                uav   = 0.5D0*(pvx(i-1,j,k)+pvx(i,j,k)) !pvx原始量
                vav   = 0.5D0*(pvy(i-1,j,k)+pvy(i,j,k))
                wav   = 0.5D0*(pvz(i-1,j,k)+pvz(i,j,k))
                q16av = 0.5D0*(q16(i-1,j,k,n)+q16(i,j,k,n))
                tav   = 0.5D0*(t(i-1,j,k)+t(i,j,k))
                call viscosity(tav,q16av,mav,kav)
                tauxx = two3*mav*(2.D0*gradfi(1,i,j,k)-gradfi(5,i,j,k)-gradfi(9,i,j,k))
                tauyy = two3*mav*(2.D0*gradfi(5,i,j,k)-gradfi(1,i,j,k)-gradfi(9,i,j,k))
                tauzz = two3*mav*(2.D0*gradfi(9,i,j,k)-gradfi(1,i,j,k)-gradfi(5,i,j,k))
                tauxy =      mav*(     gradfi(2,i,j,k)+gradfi(4,i,j,k))
                tauxz =      mav*(     gradfi(3,i,j,k)+gradfi(7,i,j,k))
                tauyz =      mav*(     gradfi(6,i,j,k)+gradfi(8,i,j,k))
                phix  = uav*tauxx + vav*tauxy + wav*tauxz  + kav*gradfi(10,i,j,k)
                phiy  = uav*tauxy + vav*tauyy + wav*tauyz  + kav*gradfi(11,i,j,k)
                phiz  = uav*tauxz + vav*tauyz + wav*tauzz  + kav*gradfi(12,i,j,k)
                flu2 = s2x(i,j,k)*tauxx + s2y(i,j,k)*tauxy+ s2z(i,j,k)*tauxz
                flu3 = s2x(i,j,k)*tauxy + s2y(i,j,k)*tauyy+ s2z(i,j,k)*tauyz
                flu4 = s2x(i,j,k)*tauxz + s2y(i,j,k)*tauyz+ s2z(i,j,k)*tauzz
                flu5 = s2x(i,j,k)*phix  + s2y(i,j,k)*phiy + s2z(i,j,k)*phiz
                qv2(i,j,k)=qv2(i,j,k)+flu2
                qv3(i,j,k)=qv3(i,j,k)+flu3
                qv4(i,j,k)=qv4(i,j,k)+flu4
                qv5(i,j,k)=qv5(i,j,k)+flu5

                qv2(i-1,j,k)=qv2(i-1,j,k)-flu2
                qv3(i-1,j,k)=qv3(i-1,j,k)-flu3
                qv4(i-1,j,k)=qv4(i-1,j,k)-flu4
                qv5(i-1,j,k)=qv5(i-1,j,k)-flu5
            end do
        end do
    end do
    !j-direction -----------------------------------------------------------------
    do k=1,nz
        do j=1,ny+1
            do i=1,nx
                uav   = 0.5D0*(pvx(i,j-1,k)+pvx(i,j,k))
                vav   = 0.5D0*(pvy(i,j-1,k)+pvy(i,j,k))
                wav   = 0.5D0*(pvz(i,j-1,k)+pvz(i,j,k))
                q16av = 0.5D0*(q16(i,j-1,k,n)+q16(i,j,k,n))
                tav   = 0.5D0*(t(i,j-1,k)+t(i,j,k))
                call viscosity(tav,q16av,mav,kav)
                tauxx = two3*mav*(2.D0*gradfj(1,i,j,k)-gradfj(5,i,j,k)-gradfj(9,i,j,k))
                tauyy = two3*mav*(2.D0*gradfj(5,i,j,k)-gradfj(1,i,j,k)-gradfj(9,i,j,k))
                tauzz = two3*mav*(2.D0*gradfj(9,i,j,k)-gradfj(1,i,j,k)-gradfj(5,i,j,k))
                tauxy =      mav*(     gradfj(2,i,j,k)+gradfj(4,i,j,k))
                tauxz =      mav*(     gradfj(3,i,j,k)+gradfj(7,i,j,k))
                tauyz =      mav*(     gradfj(6,i,j,k)+gradfj(8,i,j,k))
                phix  = uav*tauxx + vav*tauxy + wav*tauxz  + kav*gradfj(10,i,j,k)
                phiy  = uav*tauxy + vav*tauyy + wav*tauyz  + kav*gradfj(11,i,j,k)
                phiz  = uav*tauxz + vav*tauyz + wav*tauzz  + kav*gradfj(12,i,j,k)
                flu2 = s3x(i,j,k)*tauxx + s3y(i,j,k)*tauxy+ s3z(i,j,k)*tauxz
                flu3 = s3x(i,j,k)*tauxy + s3y(i,j,k)*tauyy+ s3z(i,j,k)*tauyz
                flu4 = s3x(i,j,k)*tauxz + s3y(i,j,k)*tauyz+ s3z(i,j,k)*tauzz
                flu5 = s3x(i,j,k)*phix  + s3y(i,j,k)*phiy + s3z(i,j,k)*phiz
                qv2(i,j,k)=qv2(i,j,k)+flu2
                qv3(i,j,k)=qv3(i,j,k)+flu3
                qv4(i,j,k)=qv4(i,j,k)+flu4
                qv5(i,j,k)=qv5(i,j,k)+flu5

                qv2(i,j-1,k)=qv2(i,j-1,k)-flu2
                qv3(i,j-1,k)=qv3(i,j-1,k)-flu3
                qv4(i,j-1,k)=qv4(i,j-1,k)-flu4
                qv5(i,j-1,k)=qv5(i,j-1,k)-flu5
            end do
        end do
    end do
    ! k-direction -----------------------------------------------------------------
    do k=1,nz+1
        do j=1,ny
            do i=1,nx
                uav   = 0.5D0*(pvx(i,j,k-1)+pvx(i,j,k))
                vav   = 0.5D0*(pvy(i,j,k-1)+pvy(i,j,k))
                wav   = 0.5D0*(pvz(i,j,k-1)+pvz(i,j,k))
                q16av = 0.5D0*(q16(i,j,k-1,n)+q16(i,j,k,n))
                tav   = 0.5D0*(t(i,j,k-1)+t(i,j,k))
                call viscosity(tav,q16av,mav,kav)
                tauxx = two3*mav*(2.D0*gradfk(1,i,j,k)-gradfk(5,i,j,k)-gradfk(9,i,j,k))
                tauyy = two3*mav*(2.D0*gradfk(5,i,j,k)-gradfk(1,i,j,k)-gradfk(9,i,j,k))
                tauzz = two3*mav*(2.D0*gradfk(9,i,j,k)-gradfk(1,i,j,k)-gradfk(5,i,j,k))
                tauxy =      mav*(     gradfk(2,i,j,k)+gradfk(4,i,j,k))
                tauxz =      mav*(     gradfk(3,i,j,k)+gradfk(7,i,j,k))
                tauyz =      mav*(     gradfk(6,i,j,k)+gradfk(8,i,j,k))
                phix  = uav*tauxx + vav*tauxy + wav*tauxz  + kav*gradfk(10,i,j,k)
                phiy  = uav*tauxy + vav*tauyy + wav*tauyz  + kav*gradfk(11,i,j,k)
                phiz  = uav*tauxz + vav*tauyz + wav*tauzz  + kav*gradfk(12,i,j,k)
                flu2 = s1x(i,j,k)*tauxx + s1y(i,j,k)*tauxy+ s1z(i,j,k)*tauxz
                flu3 = s1x(i,j,k)*tauxy + s1y(i,j,k)*tauyy+ s1z(i,j,k)*tauyz
                flu4 = s1x(i,j,k)*tauxz + s1y(i,j,k)*tauyz+ s1z(i,j,k)*tauzz
                flu5 = s1x(i,j,k)*phix  + s1y(i,j,k)*phiy + s1z(i,j,k)*phiz
                qv2(i,j,k)=qv2(i,j,k)+flu2
                qv3(i,j,k)=qv3(i,j,k)+flu3
                qv4(i,j,k)=qv4(i,j,k)+flu4
                qv5(i,j,k)=qv5(i,j,k)+flu5

                qv2(i,j,k-1)=qv2(i,j,k-1)-flu2
                qv3(i,j,k-1)=qv3(i,j,k-1)-flu3
                qv4(i,j,k-1)=qv4(i,j,k-1)-flu4
                qv5(i,j,k-1)=qv5(i,j,k-1)-flu5
            end do
        end do
    end do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                qv3(i,j,k)=qv3(i,j,k)+rpm*q14(i,j,k,n)*vv(i,j,k)
                qv4(i,j,k)=qv4(i,j,k)-rpm*q13(i,j,k,n)*vv(i,j,k)
            end do
        end do
    end do
    end subroutine qqqv

subroutine gradsface !界面导数值计算
    use global
    implicit none
     !交界面流动参量的导数
    gradfi=0.d0
    gradfj=0.d0
    gradfk=0.d0
    call gradsfaceI(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfi(1,1:nx+1,0:ny+1,0:nz+1))
    call gradsfaceI(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfi(2,1:nx+1,0:ny+1,0:nz+1))
    call gradsfaceI(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfi(3,1:nx+1,0:ny+1,0:nz+1))
    call gradsfaceI(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfi(4,1:nx+1,0:ny+1,0:nz+1))
    call gradsfaceI(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfi(5,1:nx+1,0:ny+1,0:nz+1))
    call gradsfaceI(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfi(6,1:nx+1,0:ny+1,0:nz+1))
    call gradsfaceI(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfi(7,1:nx+1,0:ny+1,0:nz+1))
    call gradsfaceI(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfi(8,1:nx+1,0:ny+1,0:nz+1))
    call gradsfaceI(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfi(9,1:nx+1,0:ny+1,0:nz+1))
    call gradsfaceI(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfi(10,1:nx+1,0:ny+1,0:nz+1))
    call gradsfaceI(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfi(11,1:nx+1,0:ny+1,0:nz+1))
    call gradsfaceI(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfi(12,1:nx+1,0:ny+1,0:nz+1))

    call gradsfaceJ(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfj(1,0:nx+1,1:ny+1,0:nz+1))
    call gradsfaceJ(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfj(2,0:nx+1,1:ny+1,0:nz+1))
    call gradsfaceJ(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfj(3,0:nx+1,1:ny+1,0:nz+1))
    call gradsfaceJ(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfj(4,0:nx+1,1:ny+1,0:nz+1))
    call gradsfaceJ(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfj(5,0:nx+1,1:ny+1,0:nz+1))
    call gradsfaceJ(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfj(6,0:nx+1,1:ny+1,0:nz+1))
    call gradsfaceJ(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfj(7,0:nx+1,1:ny+1,0:nz+1))
    call gradsfaceJ(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfj(8,0:nx+1,1:ny+1,0:nz+1))
    call gradsfaceJ(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfj(9,0:nx+1,1:ny+1,0:nz+1))
    call gradsfaceJ(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfj(10,0:nx+1,1:ny+1,0:nz+1))
    call gradsfaceJ(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfj(11,0:nx+1,1:ny+1,0:nz+1))
    call gradsfaceJ(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfj(12,0:nx+1,1:ny+1,0:nz+1))

    call gradsfaceK(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfk(1,0:nx+1,0:ny+1,1:nz+1))
    call gradsfaceK(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfk(2,0:nx+1,0:ny+1,1:nz+1))
    call gradsfaceK(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfk(3,0:nx+1,0:ny+1,1:nz+1))
    call gradsfaceK(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfk(4,0:nx+1,0:ny+1,1:nz+1))
    call gradsfaceK(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfk(5,0:nx+1,0:ny+1,1:nz+1))
    call gradsfaceK(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfk(6,0:nx+1,0:ny+1,1:nz+1))
    call gradsfaceK(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfk(7,0:nx+1,0:ny+1,1:nz+1))
    call gradsfaceK(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfk(8,0:nx+1,0:ny+1,1:nz+1))
    call gradsfaceK(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfk(9,0:nx+1,0:ny+1,1:nz+1))
    call gradsfaceK(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfk(10,0:nx+1,0:ny+1,1:nz+1))
    call gradsfaceK(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfk(11,0:nx+1,0:ny+1,1:nz+1))
    call gradsfaceK(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfk(12,0:nx+1,0:ny+1,1:nz+1))
    end subroutine gradsface

subroutine gradsfaceI(si,sj,sk,q,Idqd) !界面导数值计算
    use global
    implicit none
    real(8),intent(in)  ::si(0:nx+2,1:ny,1:nz),sj(1:nx,0:ny+2,1:nz),sk(1:nx,1:ny,0:nz+2),q(0:nx+1,0:ny+1,0:nz+1)
    real(8),intent(out) ::Idqd(1:nx+1,0:ny+1,0:nz+1)
    real(8) ::sx,fg,qav,rvol,sx1,sx2,qav1,qav2

    Idqd=0.d0
    do k=1,nz
        do j=1,ny
            do i=1,nx    ! 左右面通量x
                sx= 0.5D0*(si(i,j,k)+si(i+1,j,k))
                fg=q(i,j,k)*sx

                Idqd(i,j,k)  = Idqd(i,j,k)  - fg
                Idqd(i+1,j,k)= Idqd(i+1,j,k)+ fg
            end do
            sx  =  si(1,j,k)
            fg  =  0.5D0*(q(0,j,k)+q(1,j,k))*sx
            Idqd(1,j,k) = Idqd(1,j,k) + fg
            sx  =  si(nx+1,j,k)
            fg  =  0.5D0*(q(nx,j,k)+q(nx+1,j,k))*sx
            Idqd(nx+1,j,k) = Idqd(nx+1,j,k) - fg
         end do
    end do

    do k=1,nz              !上下面通量y
        do j=1,ny+1
            do i=2,nx
                sx1=0.5D0*sj(i,j,k)
                qav1=0.5D0*(q(i,j,k)+q(i,j-1,k))
                sx2=0.5D0*sj(i-1,j,k)
                qav2=0.5D0*(q(i-1,j,k)+q(i-1,j-1,k))
                fg=sx1*qav1+sx2*qav2
                Idqd(i,j,k)  = Idqd(i,j,k)  + fg
                Idqd(i,j-1,k)= Idqd(i,j-1,k)- fg
            end do
            sx  = sj(1,j,k)
            qav = 0.5D0*(q(1,j,k)+q(1,j-1,k))
            fg  = qav*sx
            Idqd(1,j,k)   = Idqd(1,j,k)   + fg
            Idqd(1,j-1,k) = Idqd(1,j-1,k) - fg
            sx  = sj(nx,j,k)
            qav = 0.5D0*(q(nx,j,k)+q(nx,j-1,k))
            fg  = qav*sx
            Idqd(nx+1,j,k)   = Idqd(nx+1,j,k)   + fg
            Idqd(nx+1,j-1,k) = Idqd(nx+1,j-1,k) - fg
        end do
    end do

    do j=1,ny
        do k=1,nz+1
            do i=2,nx
                sx1=0.5D0*sk(i,j,k)
                qav1=0.5D0*(q(i,j,k)+q(i,j,k-1))
                sx2=0.5D0*sk(i-1,j,k)
                qav2=0.5D0*(q(i-1,j,k)+q(i-1,j,k-1))
                fg=sx1*qav1+sx2*qav2
                Idqd(i,j,k)  = Idqd(i,j,k)  + fg
                Idqd(i,j,k-1)= Idqd(i,j,k-1)- fg
            end do
            sx  = sk(1,j,k)
            qav = 0.5D0*(q(1,j,k)+q(1,j,k-1))
            fg  = qav*sx
            Idqd(1,j,k)   = Idqd(1,j,k)   + fg
            Idqd(1,j,k-1) = Idqd(1,j,k-1) - fg
            sx  = sk(nx,j,k)
            qav = 0.5D0*(q(nx,j,k)+q(nx,j,k-1))
            fg  = qav*sx
            Idqd(nx+1,j,k)   = Idqd(nx+1,j,k)   + fg
            Idqd(nx+1,j,k-1) = Idqd(nx+1,j,k-1) - fg
        end do
    end do

    do k=1,nz
        do j=1,ny
            do i=2,nx
                rvol       = 2.D0/(vv(i,j,k)+vv(i-1,j,k))
                Idqd(i,j,k)= Idqd(i,j,k)*rvol
            end do
            rvol       = 1.D0/vv(1,j,k)
            Idqd(1,j,k)= Idqd(1,j,k)*rvol

            rvol          = 1.D0/vv(nx,j,k)
            Idqd(nx+1,j,k)= Idqd(nx+1,j,k)*rvol
        end do
    end do
    end subroutine gradsfaceI

subroutine gradsfaceJ(si,sj,sk,q,Jdqd) !界面导数值计算
    use global
    implicit none
    real(8),intent(in)  ::si(0:nx+2,1:ny,1:nz),sj(1:nx,0:ny+2,1:nz),sk(1:nx,1:ny,0:nz+2),q(0:nx+1,0:ny+1,0:nz+1)
    real(8),intent(out) ::Jdqd(0:nx+1,1:ny+1,0:nz+1)
    real(8) ::sx,fg,qav,rvol,sx1,sx2,qav1,qav2

    Jdqd=0.d0
    do k=1,nz
        do i=1,nx
            do j=1,ny
                sx   =  0.5D0*(sj(i,j,k)+sj(i,j+1,k))
                fg   =  q(i,j,k)*sx
                Jdqd(i,j,k)   = Jdqd(i,j,k)   - fg
                Jdqd(i,j+1,k) = Jdqd(i,j+1,k) + fg
            end do
            sx  =  sj(i,1,k)
            fg  =  0.5D0*( q(i,0,k)+ q(i,1,k))*sx
            Jdqd(i,1,k) =  Jdqd(i,1,k)  +fg
            sx  = sj(i,ny+1,k)
            fg  =  0.5D0*( q(i,ny,k)+ q(i,ny+1,k))*sx
            Jdqd(i,ny+1,k) =  Jdqd(i,ny+1,k)  -fg
         end do
    end do

    do k=1,nz
        do i=1,nx+1   !左右面通量x
            do j=2,ny
                sx1=0.5D0*si(i,j,k)
                qav1=0.5D0*(q(i,j,k)+q(i-1,j,k))
                sx2=0.5D0*si(i,j-1,k)
                qav2=0.5D0*(q(i,j-1,k)+q(i-1,j-1,k))
                fg=sx1*qav1+sx2*qav2
                Jdqd(i,j,k)   = Jdqd(i,j,k)   + fg
                Jdqd(i-1,j,k) = Jdqd(i-1,j,k) - fg
            end do
            sx  = si(i,1,k)
            qav = 0.5D0*(q(i,1,k)+q(i-1,1,k))
            fg  = qav*sx
            Jdqd(i,1,k)   = Jdqd(i,1,k)   + fg
            Jdqd(i-1,1,k) = Jdqd(i-1,1,k) - fg

            sx  = si(i,ny,k)
            qav = 0.5D0*(q(i,ny,k)+q(i-1,ny,k))
            fg  = qav*sx
            Jdqd(i,ny+1,k)   = Jdqd(i,ny+1,k)   + fg
            Jdqd(i-1,ny+1,k) = Jdqd(i-1,ny+1,k) - fg
        end do
    end do

    do i=1,nx
        do k=1,nz+1   !前后面通量z
            do j=2,ny
                sx1=0.5D0*sk(i,j,k)
                qav1=0.5D0*(q(i,j,k)+q(i,j,k-1))
                sx2=0.5D0*sk(i,j-1,k)
                qav2=0.5D0*(q(i,j-1,k)+q(i,j-1,k-1))
                fg=sx1*qav1+sx2*qav2
                Jdqd(i,j,k)   = Jdqd(i,j,k)   + fg
                Jdqd(i,j,k-1) = Jdqd(i,j,k-1) - fg
            end do
            sx  = sk(i,1,k)
            qav = 0.5D0*(q(i,1,k)+q(i,1,k-1))
            fg  = qav*sx
            Jdqd(i,1,k)   = Jdqd(i,1,k)   + fg
            Jdqd(i,1,k-1) = Jdqd(i,1,k-1) - fg

            sx  = sk(i,ny,k)
            qav = 0.5D0*(q(i,ny,k)+q(i,ny,k-1))
            fg  = qav*sx
            Jdqd(i,ny+1,k)   = Jdqd(i,ny+1,k)   + fg
            Jdqd(i,ny+1,k-1) = Jdqd(i,ny+1,k-1) - fg
        end do
    end do

    do k=1,nz
        do i=1,nx
            do j=2,ny
                rvol     = 2.D0/(vv(i,j,k)+vv(i,j-1,k))
                Jdqd(i,j,k)= Jdqd(i,j,k)*rvol
            end do
            rvol        = 1.D0/vv(i,1,k)
            Jdqd(i,1,k)   = Jdqd(i,1,k)*rvol

            rvol        = 1.D0/vv(i,ny,k)
            Jdqd(i,ny+1,k)  = Jdqd(i,ny+1,k)*rvol
        end do
    end do
   end subroutine gradsfaceJ

subroutine gradsfaceK(si,sj,sk,q,Kdqd)
    use global
    implicit none
    real(8),intent(in)  ::si(0:nx+2,1:ny,1:nz),sj(1:nx,0:ny+2,1:nz),sk(1:nx,1:ny,0:nz+2),q(0:nx+1,0:ny+1,0:nz+1)
    real(8),intent(out) ::Kdqd(0:nx+1,0:ny+1,1:nz+1)
    real(8) ::sx,fg,qav,rvol,sx1,sx2,qav1,qav2

    Kdqd=0.d0
    do j=1,ny
        do i=1,nx
            do k=1,nz
                sx= 0.5D0*(sk(i,j,k)+sk(i,j,k+1))
                fg=q(i,j,k)*sx

                Kdqd(i,j,k)  = Kdqd(i,j,k)  - fg
                Kdqd(i,j,k+1)= Kdqd(i,j,k+1)+ fg
            end do
            sx  =  sk(i,j,1)
            fg=0.5*(q(i,j,0)+q(i,j,1))*sx
            Kdqd(i,j,1) = Kdqd(i,j,1) + fg
            sx  =  sk(i,j,nz+1)
            fg=0.5*(q(i,j,nz)+q(i,j,nz+1))*sx
            Kdqd(i,j,nz+1) = Kdqd(i,j,nz+1) - fg
         end do
    end do

    do j=1,ny
        do i=1,nx+1
            do k=2,nz
                sx1=0.5D0*si(i,j,k)
                qav1=0.5D0*(q(i,j,k)+q(i-1,j,k))
                sx2=0.5D0*si(i,j,k-1)
                qav2=0.5D0*(q(i,j,k-1)+q(i-1,j,k-1))
                fg=sx1*qav1+sx2*qav2
                Kdqd(i,j,k)   = Kdqd(i,j,k)   + fg
                Kdqd(i-1,j,k) = Kdqd(i-1,j,k) - fg
            end do
            sx  = si(i,j,1)
            qav = 0.5D0*(q(i,j,1)+q(i-1,j,1))
            fg  = qav*sx
            Kdqd(i,j,1)   = Kdqd(i,j,1)   + fg
            Kdqd(i-1,j,1) = Kdqd(i-1,j,1) - fg

            sx  = si(i,j,nz)
            qav = 0.5D0*(q(i,j,nz)+q(i-1,j,nz))
            fg  = qav*sx
            Kdqd(i,j,nz+1)   = Kdqd(i,j,nz+1)   + fg
            Kdqd(i-1,j,nz+1) = Kdqd(i-1,j,nz+1) - fg
        end do
    end do

    do i=1,nx              !上下面通量y
        do j=1,ny+1
            do k=2,nz
                sx1=0.5D0*sj(i,j,k)
                qav1=0.5D0*(q(i,j,k)+q(i,j-1,k))
                sx2=0.5D0*sj(i,j,k-1)
                qav2=0.5D0*(q(i,j,k-1)+q(i,j-1,k-1))
                fg=sx1*qav1+sx2*qav2
                Kdqd(i,j,k)  = Kdqd(i,j,k)  + fg
                Kdqd(i,j-1,k)= Kdqd(i,j-1,k)- fg
            end do
            sx  = sj(i,j,1)
            qav = 0.5D0*(q(i,j,1)+q(i,j-1,1))
            fg  = qav*sx
            Kdqd(i,j,1)   = Kdqd(i,j,1)   + fg
            Kdqd(i,j-1,1) = Kdqd(i,j-1,1) - fg
            sx  = sj(i,j,nz)
            qav = 0.5D0*(q(i,j,nz)+q(i,j-1,nz))
            fg  = qav*sx
            Kdqd(i,j,nz+1)   = Kdqd(i,j,nz+1)   + fg
            Kdqd(i,j-1,nz+1) = Kdqd(i,j-1,nz+1) - fg
        end do
    end do

    do j=1,ny
        do i=1,nx
            do k=2,nz
                rvol       = 2.D0/(vv(i,j,k)+vv(i,j,k-1))
                Kdqd(i,j,k)= Kdqd(i,j,k)*rvol
            end do
            rvol       = 1.D0/vv(i,j,1)
            Kdqd(i,j,1)= Kdqd(i,j,1)*rvol

            rvol          = 1.D0/vv(i,j,nz)
            Kdqd(i,j,nz+1)= Kdqd(i,j,nz+1)*rvol
        end do
    end do
    end subroutine gradsfaceK

subroutine viscosity(temp,q6,cv,kc) !计算粘性系数和导热系数
    use global
    implicit none
    real(8),intent(in)  ::temp,q6
    real(8),intent(out) ::cv,kc
    real(8) ::cvl,cvt,fv1,tem

    cvl=cvl0*((temp/t0)**1.5)*(t0+ts)/(temp+ts)
    tem=q6/cvl
    fv1=1.d0/(1.d0+(cv1/tem)**3)
    cvt=q6*fv1
    cv=cvl+cvt
    kc=cp*(cvl/prl+cvt/prt)
    end subroutine viscosity

subroutine qqqvsa !SA方程粘性通量的计算
    use global
    implicit none
    real(8) ::flu6,q16av,tav,uq,cvl
    real(8),allocatable ::tur(:,:,:)

    allocate(tur(0:nx+1,0:ny+1,0:nz+1))
    do k=0,nz+1
        do j=0,ny+1
            do i=0,nx+1
                tur(i,j,k)=q16(i,j,k,n)/q11(i,j,k,n)
            end do
        end do
    end do
    gradfi=0.d0 !交界面流动参量的导数
    gradfj=0.d0
    gradfk=0.d0

    call gradsfaceI(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfi(13,1:nx+1,0:ny+1,0:nz+1))
    call gradsfaceI(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfi(14,1:nx+1,0:ny+1,0:nz+1))
    call gradsfaceI(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfi(15,1:nx+1,0:ny+1,0:nz+1))

    call gradsfaceJ(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfj(13,0:nx+1,1:ny+1,0:nz+1))
    call gradsfaceJ(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfj(14,0:nx+1,1:ny+1,0:nz+1))
    call gradsfaceJ(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfj(15,0:nx+1,1:ny+1,0:nz+1))

    call gradsfaceK(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfk(13,0:nx+1,0:ny+1,1:nz+1))
    call gradsfaceK(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfk(14,0:nx+1,0:ny+1,1:nz+1))
    call gradsfaceK(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfk(15,0:nx+1,0:ny+1,1:nz+1))

    qv6=0.d0
    ! i-direction -----------------------------------------------------------------
    do k=1,nz
        do j=1,ny
            do i=1,nx+1
                tav   = 0.5D0*(t(i-1,j,k)+t(i,j,k))
                cvl   = cvl0*(tav/t0)**1.5*(t0+ts)/(tav+ts)
                q16av=0.5D0*(q16(i-1,j,k,n)+q16(i,j,k,n))
                flu6  = (cvl+q16av)*(gradfi(13,i,j,k)*s2x(i,j,k)+gradfi(14,i,j,k)*s2y(i,j,k)+gradfi(15,i,j,k)*s2z(i,j,k))/sigmav
                qv6(i,j,k) =qv6(i,j,k) +flu6
                qv6(i-1,j,k)=qv6(i-1,j,k)-flu6
            end do
        end do
    end do
    ! j-direction -----------------------------------------------------------------
    do k=1,nz
        do j=1,ny+1
            do i=1,nx
                tav   = 0.5D0*(t(i,j-1,k)+t(i,j,k))
                cvl   = cvl0*(tav/t0)**1.5*(t0+ts)/(tav+ts)
                q16av=0.5D0*(q16(i,j-1,k,n)+q16(i,j,k,n))
                flu6  = (cvl+q16av)*(gradfj(13,i,j,k)*s3x(i,j,k)+gradfj(14,i,j,k)*s3y(i,j,k)+gradfj(15,i,j,k)*s3z(i,j,k))/sigmav
                qv6(i,j,k) =qv6(i,j,k) +flu6
                qv6(i,j-1,k)=qv6(i,j-1,k)-flu6
            end do
        end do
    end do
    ! k-direction -----------------------------------------------------------------
    do k=1,nz+1
        do j=1,ny
            do i=1,nx
                tav   = 0.5D0*(t(i,j,k-1)+t(i,j,k))
                cvl   = cvl0*(tav/t0)**1.5*(t0+ts)/(tav+ts)
                q16av=0.5D0*(q16(i,j,k-1,n)+q16(i,j,k,n))
                flu6  = (cvl+q16av)*(gradfk(13,i,j,k)*s1x(i,j,k)+gradfk(14,i,j,k)*s1y(i,j,k)+gradfk(15,i,j,k)*s1z(i,j,k))/sigmav
                qv6(i,j,k) =qv6(i,j,k) +flu6
                qv6(i,j,k-1)=qv6(i,j,k-1)-flu6
            end do
        end do
    end do
    deallocate(tur)
    call SAsource
    end subroutine qqqvsa

subroutine SAsource !求网格中心的方程六的源项
    use global
    implicit none
    real(8),allocatable ::tur(:,:,:),pwx(:,:,:),pwy(:,:,:),pwz(:,:,:)
    real(8) ::cvl,tem,fv1,fv2,fv3,rp1,rp2,rp3,w12,w13,w23,w21,w31,w32,ww2,ww,svot,vm,gv
    real(8) ::s11,s22,s33,s12,s13,s23,ss2,ss,dd,fr1r1,fr1r2,fr1,yv,ra,ga,fw,dv

    allocate(tur(0:nx+1,0:ny+1,0:nz+1))
    allocate(pwx(0:nx+1,0:ny+1,0:nz+1))
    allocate(pwy(0:nx+1,0:ny+1,0:nz+1))
    allocate(pwz(0:nx+1,0:ny+1,0:nz+1))

    yy0(0,1:ny,1:nz)=yy0(1,1:ny,1:nz) !网格中心点坐标，虚拟网格的处理
    zz0(0,1:ny,1:nz)=zz0(1,1:ny,1:nz)
    yy0(nx+1,1:ny,1:nz)=yy0(nx,1:ny,1:nz)
    zz0(nx+1,1:ny,1:nz)=zz0(nx,1:ny,1:nz)

    yy0(1:nx,0,1:nz)=yy0(1:nx,1,1:nz)
    zz0(1:nx,0,1:nz)=zz0(1:nx,1,1:nz)
    yy0(1:nx,ny+1,1:nz)=yy0(1:nx,ny,1:nz)
    zz0(1:nx,ny+1,1:nz)=zz0(1:nx,ny,1:nz)

    yy0(1:nx,1:ny,0)=yy0(1:nx,1:ny,1)
    zz0(1:nx,1:ny,0)=zz0(1:nx,1:ny,1)
    yy0(1:nx,1:ny,nz+1)=yy0(1:nx,1:ny,nz)
    zz0(1:nx,1:ny,nz+1)=zz0(1:nx,1:ny,nz)
    do k=0,nz+1
        do j=0,ny+1
            do i=0,nx+1
                tur(i,j,k)=q16(i,j,k,n)/q11(i,j,k,n)
                pwx(i,j,k)=pvx(i,j,k)
                pwy(i,j,k)=pvy(i,j,k)+rpm*zz0(i,j,k)
                pwz(i,j,k)=pvz(i,j,k)-rpm*yy0(i,j,k)
            end do
        end do
    end do
    !********网格中心导数***********
    gradc=0.d0 !交界面流动参量的导数
    call gradscentre(1,pwx(0:nx+1,0:ny+1,0:nz+1),gradc(1,0:nx+1,0:ny+1,0:nz+1))
    call gradscentre(2,pwx(0:nx+1,0:ny+1,0:nz+1),gradc(2,0:nx+1,0:ny+1,0:nz+1))
    call gradscentre(3,pwx(0:nx+1,0:ny+1,0:nz+1),gradc(3,0:nx+1,0:ny+1,0:nz+1))
    call gradscentre(1,pwy(0:nx+1,0:ny+1,0:nz+1),gradc(4,0:nx+1,0:ny+1,0:nz+1))
    call gradscentre(2,pwy(0:nx+1,0:ny+1,0:nz+1),gradc(5,0:nx+1,0:ny+1,0:nz+1))
    call gradscentre(3,pwy(0:nx+1,0:ny+1,0:nz+1),gradc(6,0:nx+1,0:ny+1,0:nz+1))
    call gradscentre(1,pwz(0:nx+1,0:ny+1,0:nz+1),gradc(7,0:nx+1,0:ny+1,0:nz+1))
    call gradscentre(2,pwz(0:nx+1,0:ny+1,0:nz+1),gradc(8,0:nx+1,0:ny+1,0:nz+1))
    call gradscentre(3,pwz(0:nx+1,0:ny+1,0:nz+1),gradc(9,0:nx+1,0:ny+1,0:nz+1))
    call gradscentre(1,tur(0:nx+1,0:ny+1,0:nz+1),gradc(10,0:nx+1,0:ny+1,0:nz+1))
    call gradscentre(2,tur(0:nx+1,0:ny+1,0:nz+1),gradc(11,0:nx+1,0:ny+1,0:nz+1))
    call gradscentre(3,tur(0:nx+1,0:ny+1,0:nz+1),gradc(12,0:nx+1,0:ny+1,0:nz+1))
    gradcs=0.d0  !交界面流动参量的导数
    call dsdt(gradc(1,0:nx+1,0:ny+1,0:nz+1),gradcs(1,0:nx+1,0:ny+1,0:nz+1))
    call dsdt(gradc(2,0:nx+1,0:ny+1,0:nz+1),gradcs(2,0:nx+1,0:ny+1,0:nz+1))
    call dsdt(gradc(3,0:nx+1,0:ny+1,0:nz+1),gradcs(3,0:nx+1,0:ny+1,0:nz+1))
    call dsdt(gradc(4,0:nx+1,0:ny+1,0:nz+1),gradcs(4,0:nx+1,0:ny+1,0:nz+1))
    call dsdt(gradc(5,0:nx+1,0:ny+1,0:nz+1),gradcs(5,0:nx+1,0:ny+1,0:nz+1))
    call dsdt(gradc(6,0:nx+1,0:ny+1,0:nz+1),gradcs(6,0:nx+1,0:ny+1,0:nz+1))
    call dsdt(gradc(7,0:nx+1,0:ny+1,0:nz+1),gradcs(7,0:nx+1,0:ny+1,0:nz+1))
    call dsdt(gradc(8,0:nx+1,0:ny+1,0:nz+1),gradcs(8,0:nx+1,0:ny+1,0:nz+1))
    call dsdt(gradc(9,0:nx+1,0:ny+1,0:nz+1),gradcs(9,0:nx+1,0:ny+1,0:nz+1))
    do k=1,nz
        do j=1,ny
            do i=1,nx
                dv=cb2/sigmav*(gradc(10,i,j,k)**2+gradc(11,i,j,k)**2+gradc(12,i,j,k)**2)*q11(i,j,k,n)
                cvl=cvl0*(t(i,j,k)/t0)**1.5*(t0+ts)/(t(i,j,k)+ts)
                tem=q16(i,j,k,n)/cvl
                tem=max(tem,1.d-4)
                fv1=1.d0/(1.d0+(cv1/tem)**3)
                fv2=(1.d0+tem/cv2)**(-3)
                fv3=(1.d0+tem*fv1)*(1.d0-fv2)/tem
                w12=0.5d0*(gradc(2,i,j,k)-gradc(4,i,j,k))
                w13=0.5d0*(gradc(3,i,j,k)-gradc(7,i,j,k))
                w23=0.5d0*(gradc(6,i,j,k)-gradc(8,i,j,k))
                ww2=4.d0*(w12*w12+w13*w13+w23*w23)
                ww=sqrt(ww2)
                vm=tur(i,j,k)/(kap*kap*dmini(i,j,k)*dmini(i,j,k))
                svot=ww*fv3+vm*fv2
                gv=cb1*q16(i,j,k,n)*svot
                s11=gradc(1,i,j,k)
                s22=gradc(5,i,j,k)
                s33=gradc(9,i,j,k)
                s12=0.5d0*(gradc(2,i,j,k)+gradc(4,i,j,k))
                s13=0.5d0*(gradc(3,i,j,k)+gradc(7,i,j,k))
                s23=0.5d0*(gradc(6,i,j,k)+gradc(8,i,j,k))
                ss2=4.d0*(s12*s12+s13*s13+s23*s23)+2.d0*(s11*s11+s22*s22+s33*s33)
                ss=sqrt(ss2)
                rp1=rpm
                rp2=0.d0
                rp3=0.d0
                w12=0.5d0*(gradc(2,i,j,k)-gradc(4,i,j,k))-rp3
                w13=0.5d0*(gradc(3,i,j,k)-gradc(7,i,j,k))+rp2
                w23=0.5d0*(gradc(6,i,j,k)-gradc(8,i,j,k))-rp1
                w21=-w12
                w31=-w13
                w32=-w23
                ww2=4.d0*(w12*w12+w13*w13+w23*w23)
                ww=sqrt(ww2)
                fr1r1=ss/ww
                dd=max(ss2,0.09d0*tur(i,j,k)*tur(i,j,k))
                fr1r2=2.d0/dd/sqrt(dd)/ww*(gradcs(1,i,j,k)*(w12*s12+w13*s13)+(0.5d0*(gradcs(2,i,j,k)+gradcs(4,i,j,k))-s13*rp1)*(w12*s22+w13*s23+w21*s11+w23*s13)&
                      +(0.5d0*(gradcs(3,i,j,k)+gradcs(7,i,j,k))+s12*rp1)*(w12*s23+w13*s33+w31*s11+w32*s12)+(gradcs(5,i,j,k)-2.d0*s23*rp1)*(w21*s12+w23*s23)&
                      +(0.5d0*(gradcs(6,i,j,k)+gradcs(8,i,j,k))+(s22-s33)*rp1)*(w21*s13+w23*s33+w31*s12+w32*s22)+(gradcs(9,i,j,k)+2.d0*s23*rp1)*(w31*s13+w32*s23))
                fr1=(1d0+cr1)*2.d0*fr1r1/(1.d0+fr1r1)*(1.d0-cr3*atan(cr2*fr1r2))-cr1
                gv=gv*fr1
                ra=vm/svot
                ga=ra+cw2*(ra**6-ra)
                fw=((ga**(-6)+cw3**(-6))/(1.d0+cw3**(-6)))**(-1.d0/6.d0)
                yv=cw1*fw*q16(i,j,k,n)*vm*kap*kap
                qv6(i,j,k)=qv6(i,j,k)+(gv-yv+dv)*vv(i,j,k)
            end do
        end do
    end do
    deallocate(tur)
    deallocate(pwx)
    deallocate(pwy)
    deallocate(pwz)
    end subroutine SAsource

subroutine gradscentre(direction,q,dqd) !网格中心梯度
    use global
    integer,intent(in)  ::direction
    real(8),intent(in)  ::q(0:nx+1,0:ny+1,0:nz+1)
    real(8),intent(out) ::dqd(0:nx+1,0:ny+1,0:nz+1)
    real(8) ::si,flu

    dqd=0.d0
    !*********x方向
    do k=1,nz
        do j=1,ny
            do i=1,nx+1
                if (direction==1)then
                    si=s2x(i,j,k)
                else if (direction==2)then
                    si=s2y(i,j,k)
                else if (direction==3)then
                    si=s2z(i,j,k)
                end if
                flu=0.5d0*(q(i,j,k)+q(i-1,j,k))*si
                dqd(i,j,k)=dqd(i,j,k)+flu
                dqd(i-1,j,k)=dqd(i-1,j,k)-flu
            end do
        end do
    end do
  !*********y方向
    do k=1,nz
        do j=1,ny+1
            do i=1,nx
                if (direction==1)then
                    si=s3x(i,j,k)
                else if (direction==2)then
                    si=s3y(i,j,k)
                else if (direction==3)then
                    si=s3z(i,j,k)
                end if
                flu=0.5d0*(q(i,j-1,k)+q(i,j,k))*si
                dqd(i,j,k)=dqd(i,j,k)+flu
                dqd(i,j-1,k)=dqd(i,j-1,k)-flu
            end do
        end do
    end do
  !*********z方向
    do k=1,nz+1
        do j=1,ny
            do i=1,nx
                if (direction==1)then
                    si=s1x(i,j,k)
                else if (direction==2)then
                    si=s1y(i,j,k)
                else if (direction==3)then
                    si=s1z(i,j,k)
                end if
                flu=0.5d0*(q(i,j,k)+q(i,j,k-1))*si
                dqd(i,j,k)=dqd(i,j,k)+flu
                dqd(i,j,k-1)=dqd(i,j,k-1)-flu
            end do
        end do
    end do

    do k=1,nz
        do j=1,ny
            do i=1,nx
                dqd(i,j,k)=dqd(i,j,k)/vv(i,j,k)
            end do
        end do
    end do
    end subroutine gradscentre

subroutine dsdt(q,s) !旋转因子中的参数
    use global
    implicit none
    real(8),intent(in)  ::q(0:nx+1,0:ny+1,0:nz+1)
    real(8),intent(out) ::s(0:nx+1,0:ny+1,0:nz+1)
    real(8) ::vf,rf,qq1,flu

    s=0.d0
    !******************x方向********************
    do k=1,nz
        do j=1,ny
            do i=1,nx+1
               vx=0.5d0*(pvx(i,j,k)+pvx(i-1,j,k))
               vy=0.5d0*(pvy(i,j,k)+pvy(i-1,j,k))
               vz=0.5d0*(pvz(i,j,k)+pvz(i-1,j,k))
               vf=vx*s2x(i,j,k)+vy*s2y(i,j,k)+vz*s2z(i,j,k)
               rf=rpm*zz02(i,j,k)*s2y(i,j,k)-rpm*yy02(i,j,k)*s2z(i,j,k)
               vf=vf+rf
               qq1=0.5d0*(q(i,j,k)+q(i-1,j,k))
               flu=qq1*vf
               s(i,j,k)=s(i,j,k)+flu
               s(i-1,j,k)=s(i-1,j,k)-flu
            end do
        end do
    end do
   !******************y方向********************
    do k=1,nz
        do j=1,ny+1
            do i=1,nx
               vx=0.5d0*(pvx(i,j,k)+pvx(i,j-1,k))
               vy=0.5d0*(pvy(i,j,k)+pvy(i,j-1,k))
               vz=0.5d0*(pvz(i,j,k)+pvz(i,j-1,k))
               vf=vx*s3x(i,j,k)+vy*s3y(i,j,k)+vz*s3z(i,j,k)
               rf=rpm*zz03(i,j,k)*s3y(i,j,k)-rpm*yy03(i,j,k)*s3z(i,j,k)
               vf=vf+rf
               qq1=0.5d0*(q(i,j,k)+q(i,j-1,k))
               flu=qq1*vf
               s(i,j,k)=s(i,j,k)+flu
               s(i,j-1,k)=s(i,j-1,k)-flu
            end do
        end do
    end do
 !******************z方向********************
    do k=1,nz+1
        do j=1,ny
            do i=1,nx
                vx=0.5d0*(pvx(i,j,k)+pvx(i,j,k-1))
                vy=0.5d0*(pvy(i,j,k)+pvy(i,j,k-1))
                vz=0.5d0*(pvz(i,j,k)+pvz(i,j,k-1))
                vf=vx*s1x(i,j,k)+vy*s1y(i,j,k)+vz*s1z(i,j,k)
                rf=rpm*zz01(i,j,k)*s1y(i,j,k)-rpm*yy01(i,j,k)*s1z(i,j,k)
                vf=vf+rf
                qq1=0.5d0*(q(i,j,k)+q(i,j,k-1))
                flu=qq1*vf
                s(i,j,k)=s(i,j,k)+flu
                s(i,j,k-1)=s(i,j,k-1)-flu
            end do
        end do
    end do

    do k=1,nz
        do j=1,ny
            do i=1,nx
                s(i,j,k)=s(i,j,k)/vv(i,j,k)
            end do
        end do
    end do
    end subroutine dsdt

subroutine rrr   !计算方程残差
    use global
    implicit none

    rr1=0.d0
    rr2=0.d0
    rr3=0.d0
    rr4=0.d0
    rr5=0.d0
    do k=1,nz
        do j=1,ny
            do i=1,nx !qp驱动源项P2h，多重网格粗网格迭代中,qc对流项
                rr1(i,j,k)=qp1(i,j,k)-qc1(i,j,k)+av1(i,j,k)-ts1(i,j,k)*vv(i,j,k)
                rr2(i,j,k)=qp2(i,j,k)-qc2(i,j,k)+av2(i,j,k)-ts2(i,j,k)*vv(i,j,k)+qv2(i,j,k)
                rr3(i,j,k)=qp3(i,j,k)-qc3(i,j,k)+av3(i,j,k)-ts3(i,j,k)*vv(i,j,k)+qv3(i,j,k)
                rr4(i,j,k)=qp4(i,j,k)-qc4(i,j,k)+av4(i,j,k)-ts4(i,j,k)*vv(i,j,k)+qv4(i,j,k)
                rr5(i,j,k)=qp5(i,j,k)-qc5(i,j,k)+av5(i,j,k)-ts5(i,j,k)*vv(i,j,k)+qv5(i,j,k)
            end do
        end do
    end do
    end subroutine rrr

subroutine pred(ims) !每步R-K推进方法
    use global
    implicit none
    integer ::ims
    real(8) ::dtime

    call rrr
    py1=0.d0
    py2=0.d0
    py3=0.d0
    py4=0.d0
    py5=0.d0
    do k=1,nz
        do j=1,ny
            do i=1,nx
                dtime=timl/time(i,j,k)
                py1(i,j,k)=dtime*rr1(i,j,k)
                py2(i,j,k)=dtime*rr2(i,j,k)
                py3(i,j,k)=dtime*rr3(i,j,k)
                py4(i,j,k)=dtime*rr4(i,j,k)
                py5(i,j,k)=dtime*rr5(i,j,k)
            end do
        end do
    end do
    if(ims==1)then
        call ave
    end if
    do k=1,nz
        do j=1,ny
            do i=1,nx
                q11(i,j,k,n)=q01(i,j,k)+py1(i,j,k)
                q12(i,j,k,n)=q02(i,j,k)+py2(i,j,k)
                q13(i,j,k,n)=q03(i,j,k)+py3(i,j,k)
                q14(i,j,k,n)=q04(i,j,k)+py4(i,j,k)
                q15(i,j,k,n)=q05(i,j,k)+py5(i,j,k)
            end do
        end do
    end do
    end subroutine pred

subroutine predsa(ims) !SA方程每步R-K推进方法
    use global
    implicit none
    integer ::ims
    real(8) ::dtime

    py6=0.d0 !时间谱中
    do k=1,nz
        do j=1,ny
            do i=1,nx
                dtime=timl/time(i,j,k)
                py6(i,j,k)=dtime*(-qc6(i,j,k)+av6(i,j,k)-ts6(i,j,k)*vv(i,j,k)+qv6(i,j,k))
            end do
        end do
    end do
    if(ims==1)then
        call avesa
    end if
    do k=1,nz
        do j=1,ny
            do i=1,nx
                q16(i,j,k,n)=q06(i,j,k)+py6(i,j,k)
            end do
        end do
    end do
    end subroutine predsa

subroutine ave  !隐式残差光顺
    use global
    implicit none
    real(8),allocatable ::ax(:),bx(:),ay(:),by(:),az(:),bz(:)

    allocate(ax(nx))
    allocate(bx(nx))
    allocate(ay(ny))
    allocate(by(ny))
    allocate(az(nz))
    allocate(bz(nz))
    ax=-ta
    ay=-ta
    az=-ta
    bx=1.d0+2.d0*ta
    by=1.d0+2.d0*ta
    bz=1.d0+2.d0*ta
    do k=1,nz
        do j=1,ny
            call tdma(ax(2:nx),bx(1:nx),ax(1:nx-1),py1(1:nx,j,k),py1(1:nx,j,k),nx)
            call tdma(ax(2:nx),bx(1:nx),ax(1:nx-1),py2(1:nx,j,k),py2(1:nx,j,k),nx)
            call tdma(ax(2:nx),bx(1:nx),ax(1:nx-1),py3(1:nx,j,k),py3(1:nx,j,k),nx)
            call tdma(ax(2:nx),bx(1:nx),ax(1:nx-1),py4(1:nx,j,k),py4(1:nx,j,k),nx)
            call tdma(ax(2:nx),bx(1:nx),ax(1:nx-1),py5(1:nx,j,k),py5(1:nx,j,k),nx)
        end do
    end do
    do k=1,nz
        do i=1,nx
            call tdma(ay(2:ny),by(1:ny),ay(1:ny-1),py1(i,1:ny,k),py1(i,1:ny,k),ny)
            call tdma(ay(2:ny),by(1:ny),ay(1:ny-1),py2(i,1:ny,k),py2(i,1:ny,k),ny)
            call tdma(ay(2:ny),by(1:ny),ay(1:ny-1),py3(i,1:ny,k),py3(i,1:ny,k),ny)
            call tdma(ay(2:ny),by(1:ny),ay(1:ny-1),py4(i,1:ny,k),py4(i,1:ny,k),ny)
            call tdma(ay(2:ny),by(1:ny),ay(1:ny-1),py5(i,1:ny,k),py5(i,1:ny,k),ny)
        end do
    end do
    do j=1,ny
        do i=1,nx
            call tdma(az(2:nz),bz(1:nz),az(1:nz-1),py1(i,j,1:nz),py1(i,j,1:nz),nz)
            call tdma(az(2:nz),bz(1:nz),az(1:nz-1),py2(i,j,1:nz),py2(i,j,1:nz),nz)
            call tdma(az(2:nz),bz(1:nz),az(1:nz-1),py3(i,j,1:nz),py3(i,j,1:nz),nz)
            call tdma(az(2:nz),bz(1:nz),az(1:nz-1),py4(i,j,1:nz),py4(i,j,1:nz),nz)
            call tdma(az(2:nz),bz(1:nz),az(1:nz-1),py5(i,j,1:nz),py5(i,j,1:nz),nz)
        end do
    end do
    deallocate(ax)
    deallocate(bx)
    deallocate(ay)
    deallocate(by)
    deallocate(az)
    deallocate(bz)
    end subroutine ave

subroutine avesa !SA方程上的隐式残差光顺
    use global
    implicit none
    real(8),allocatable ::ax(:),bx(:),ay(:),by(:),az(:),bz(:)

    allocate(ax(nx))
    allocate(bx(nx))
    allocate(ay(ny))
    allocate(by(ny))
    allocate(az(nz))
    allocate(bz(nz))
    ax=-ta
    ay=-ta
    az=-ta
    bx=1.d0+2.d0*ta
    by=1.d0+2.d0*ta
    bz=1.d0+2.d0*ta
    do k=1,nz
        do j=1,ny
            call tdma(ax(2:nx),bx(1:nx),ax(1:nx-1),py6(1:nx,j,k),py6(1:nx,j,k),nx)
        end do
    end do
    do k=1,nz
        do i=1,nx
            call tdma(ay(2:ny),by(1:ny),ay(1:ny-1),py6(i,1:ny,k),py6(i,1:ny,k),ny)
        end do
    end do
    do j=1,ny
        do i=1,nx
            call tdma(az(2:nz),bz(1:nz),az(1:nz-1),py6(i,j,1:nz),py6(i,j,1:nz),nz)
        end do
    end do
    deallocate(ax)
    deallocate(bx)
    deallocate(ay)
    deallocate(by)
    deallocate(az)
    deallocate(bz)
    end subroutine avesa

subroutine tdma(a,b,c,d,x,n) !TDMA算法,tdma算法求解三对角矩阵
    implicit none
    integer,intent(in) ::n
    real(8),intent(in) ::a(2:n),b(1:n),c(1:n-1),d(1:n)
    real(8),intent(out)::x(1:n)
    integer ::i
    real(8) ::u(1:n),l(2:n),y(1:n)

    u(1)=b(1)
    y(1)=d(1)
    do i=2,n
        l(i)=a(i)/u(i-1)
        u(i)=b(i)-l(i)*c(i-1)
        y(i)=d(i)-l(i)*y(i-1)
    end do
    x(n)=y(n)/u(n)
    do i=n-1,1,-1
        x(i)=(y(i)-c(i)*x(i+1))/u(i)
    end do
end subroutine tdma

subroutine residual !迭代中最大残差量的计算
    use global
    implicit none

    rms(n)=0.d0
    do k=1,nz
        do j=1,ny
            do i=1,nx
                rms(n)=rms(n)+((Q11(i,j,k,n)-Q01(i,j,k))*time(i,j,k)/vv(i,j,k)/cfl)**2
            end do
        end do
    end do
    rms(n)=0.5d0*log10(rms(n)/dble(nx*ny*nz))
    if (rms(n)>rmsm)then
        rmsm=rms(n)
    end if
end subroutine residual

subroutine init !限制,多重网格中的限制算子
    use global
    implicit none
    integer ::i1,j1,k1
    real(8) ::val1,val2,val3,val4,val5,val6,val7,val8,vall

    do k=1,nz
        k1=2*k !每层网格有重新赋值，解和残差限制到粗网格上
        do j=1,ny
            j1=2*j
            do i=1,nx
                i1=2*i
                val1=vvn(i1,j1,k1,ign+1)
                val2=vvn(i1-1,j1,k1,ign+1)
                val3=vvn(i1,j1,k1-1,ign+1)
                val4=vvn(i1-1,j1,k1-1,ign+1)
                val5=vvn(i1,j1-1,k1,ign+1)
                val6=vvn(i1-1,j1-1,k1,ign+1)
                val7=vvn(i1,j1-1,k1-1,ign+1)
                val8=vvn(i1-1,j1-1,k1-1,ign+1)
                vall=val1+val2+val3+val4+val5+val6+val7+val8
                q11(i,j,k,:)=(q11(i1,j1,k1,:)*val1+q11(i1-1,j1,k1,:)*val2+q11(i1,j1,k1-1,:)*val3+q11(i1-1,j1,k1-1,:)*val4+                &
                              q11(i1,j1-1,k1,:)*val5+q11(i1-1,j1-1,k1,:)*val6+q11(i1,j1-1,k1-1,:)*val7+q11(i1-1,j1-1,k1-1,:)*val8)/vall
                q12(i,j,k,:)=(q12(i1,j1,k1,:)*val1+q12(i1-1,j1,k1,:)*val2+q12(i1,j1,k1-1,:)*val3+q12(i1-1,j1,k1-1,:)*val4+                &
                              q12(i1,j1-1,k1,:)*val5+q12(i1-1,j1-1,k1,:)*val6+q12(i1,j1-1,k1-1,:)*val7+q12(i1-1,j1-1,k1-1,:)*val8)/vall
                q13(i,j,k,:)=(q13(i1,j1,k1,:)*val1+q13(i1-1,j1,k1,:)*val2+q13(i1,j1,k1-1,:)*val3+q13(i1-1,j1,k1-1,:)*val4+                &
                              q13(i1,j1-1,k1,:)*val5+q13(i1-1,j1-1,k1,:)*val6+q13(i1,j1-1,k1-1,:)*val7+q13(i1-1,j1-1,k1-1,:)*val8)/vall
                q14(i,j,k,:)=(q14(i1,j1,k1,:)*val1+q14(i1-1,j1,k1,:)*val2+q14(i1,j1,k1-1,:)*val3+q14(i1-1,j1,k1-1,:)*val4+                &
                              q14(i1,j1-1,k1,:)*val5+q14(i1-1,j1-1,k1,:)*val6+q14(i1,j1-1,k1-1,:)*val7+q14(i1-1,j1-1,k1-1,:)*val8)/vall
                q15(i,j,k,:)=(q15(i1,j1,k1,:)*val1+q15(i1-1,j1,k1,:)*val2+q15(i1,j1,k1-1,:)*val3+q15(i1-1,j1,k1-1,:)*val4+                &
                              q15(i1,j1-1,k1,:)*val5+q15(i1-1,j1-1,k1,:)*val6+q15(i1,j1-1,k1-1,:)*val7+q15(i1-1,j1-1,k1-1,:)*val8)/vall
                q16(i,j,k,:)=(q16(i1,j1,k1,:)*val1+q16(i1-1,j1,k1,:)*val2+q16(i1,j1,k1-1,:)*val3+q16(i1-1,j1,k1-1,:)*val4+                &
                              q16(i1,j1-1,k1,:)*val5+q16(i1-1,j1-1,k1,:)*val6+q16(i1,j1-1,k1-1,:)*val7+q16(i1-1,j1-1,k1-1,:)*val8)/vall
                rr1(i,j,k)=rr1(i1,j1,k1)+rr1(i1-1,j1,k1)+rr1(i1,j1,k1-1)+rr1(i1-1,j1,k1-1)+                &
                           rr1(i1,j1-1,k1)+rr1(i1-1,j1-1,k1)+rr1(i1,j1-1,k1-1)+rr1(i1-1,j1-1,k1-1)
                rr2(i,j,k)=rr2(i1,j1,k1)+rr2(i1-1,j1,k1)+rr2(i1,j1,k1-1)+rr2(i1-1,j1,k1-1)+                &
                           rr2(i1,j1-1,k1)+rr2(i1-1,j1-1,k1)+rr2(i1,j1-1,k1-1)+rr2(i1-1,j1-1,k1-1)
                rr3(i,j,k)=rr3(i1,j1,k1)+rr3(i1-1,j1,k1)+rr3(i1,j1,k1-1)+rr3(i1-1,j1,k1-1)+                &
                           rr3(i1,j1-1,k1)+rr3(i1-1,j1-1,k1)+rr3(i1,j1-1,k1-1)+rr3(i1-1,j1-1,k1-1)
                rr4(i,j,k)=rr4(i1,j1,k1)+rr4(i1-1,j1,k1)+rr4(i1,j1,k1-1)+rr4(i1-1,j1,k1-1)+                &
                           rr4(i1,j1-1,k1)+rr4(i1-1,j1-1,k1)+rr4(i1,j1-1,k1-1)+rr4(i1-1,j1-1,k1-1)
                rr5(i,j,k)=rr5(i1,j1,k1)+rr5(i1-1,j1,k1)+rr5(i1,j1,k1-1)+rr5(i1-1,j1,k1-1)+                &
                           rr5(i1,j1-1,k1)+rr5(i1-1,j1-1,k1)+rr5(i1,j1-1,k1-1)+rr5(i1-1,j1-1,k1-1)
            end do
        end do
     end do
end subroutine init

subroutine copr !粗网格驱动项P2h
    use global
    implicit none

    call tsd
    call ppp
    call bc
    call step
    call dddc
    call qqq
    call qqqv

    do k=1,nz
        do j=1,ny
            do i=1,nx
                qp1(i,j,k)=rr1(i,j,k)+qc1(i,j,k)-av1(i,j,k)+ts1(i,j,k)*vv(i,j,k)
                qp2(i,j,k)=rr2(i,j,k)+qc2(i,j,k)-av2(i,j,k)+ts2(i,j,k)*vv(i,j,k)-qv2(i,j,k)
                qp3(i,j,k)=rr3(i,j,k)+qc3(i,j,k)-av3(i,j,k)+ts3(i,j,k)*vv(i,j,k)-qv3(i,j,k)
                qp4(i,j,k)=rr4(i,j,k)+qc4(i,j,k)-av4(i,j,k)+ts4(i,j,k)*vv(i,j,k)-qv4(i,j,k)
                qp5(i,j,k)=rr5(i,j,k)+qc5(i,j,k)-av5(i,j,k)+ts5(i,j,k)*vv(i,j,k)-qv5(i,j,k)
            end do
        end do
    end do
    end subroutine copr

subroutine update !插值
    use global
    implicit none
    integer ::ig,i1,j1,k1,iw,ie,jn,js,kf,kb
    real,external ::  dp
    real(8) :: dpy,dpiw,dpie,dpjn,dpjs,dpkf,dpkb
    do ig=ign,nng-1
        nx=nnx(ig)
        ny=nny(ig)
        nz=nnz(ig)
        do k=1,nz
            k1=2*k-1
            kf=max(1,k-1)
            kb=min(nz,k+1)
            do j=1,ny
                j1=2*j-1
                js=max(1,j-1)
                jn=min(ny,j+1)
                do i=1,nx
                    i1=2*i-1
                    iw=max(1,i-1)
                    ie=min(nx,i+1)

                    dpy=py1(i,j,k)
                    dpiw=dp(py1(iw,j,k),dpy,py1(ie,j,k))
                    dpie=dp(py1(ie,j,k),dpy,py1(iw,j,k))
                    dpjs=dp(py1(i,js,k),dpy,py1(i,jn,k))
                    dpjn=dp(py1(i,jn,k),dpy,py1(i,js,k))
                    dpkf=dp(py1(i,j,kf),dpy,py1(i,j,kb))
                    dpkb=dp(py1(i,j,kb),dpy,py1(i,j,kf))
                    q01(i1,j1,k1)=dpy+dpiw+dpjs+dpkf
                    q01(i1+1,j1,k1)=dpy+dpie+dpjs+dpkf
                    q01(i1,j1,k1+1)=dpy+dpiw+dpjs+dpkb
                    q01(i1+1,j1,k1+1)=dpy+dpie+dpjs+dpkb
                    q01(i1,j1+1,k1)=dpy+dpiw+dpjn+dpkf
                    q01(i1+1,j1+1,k1)=dpy+dpie+dpjn+dpkf
                    q01(i1,j1+1,k1+1)=dpy+dpiw+dpjn+dpkb
                    q01(i1+1,j1+1,k1+1)=dpy+dpie+dpjn+dpkb

                    dpy=py2(i,j,k)
                    dpiw=dp(py2(iw,j,k),dpy,py2(ie,j,k))
                    dpie=dp(py2(ie,j,k),dpy,py2(iw,j,k))
                    dpjs=dp(py2(i,js,k),dpy,py2(i,jn,k))
                    dpjn=dp(py2(i,jn,k),dpy,py2(i,js,k))
                    dpkf=dp(py2(i,j,kf),dpy,py2(i,j,kb))
                    dpkb=dp(py2(i,j,kb),dpy,py2(i,j,kf))
                    q02(i1,j1,k1)=dpy+dpiw+dpjs+dpkf
                    q02(i1+1,j1,k1)=dpy+dpie+dpjs+dpkf
                    q02(i1,j1,k1+1)=dpy+dpiw+dpjs+dpkb
                    q02(i1+1,j1,k1+1)=dpy+dpie+dpjs+dpkb
                    q02(i1,j1+1,k1)=dpy+dpiw+dpjn+dpkf
                    q02(i1+1,j1+1,k1)=dpy+dpie+dpjn+dpkf
                    q02(i1,j1+1,k1+1)=dpy+dpiw+dpjn+dpkb
                    q02(i1+1,j1+1,k1+1)=dpy+dpie+dpjn+dpkb

                    dpy=py3(i,j,k)
                    dpiw=dp(py3(iw,j,k),dpy,py3(ie,j,k))
                    dpie=dp(py3(ie,j,k),dpy,py3(iw,j,k))
                    dpjs=dp(py3(i,js,k),dpy,py3(i,jn,k))
                    dpjn=dp(py3(i,jn,k),dpy,py3(i,js,k))
                    dpkf=dp(py3(i,j,kf),dpy,py3(i,j,kb))
                    dpkb=dp(py3(i,j,kb),dpy,py3(i,j,kf))
                    q03(i1,j1,k1)=dpy+dpiw+dpjs+dpkf
                    q03(i1+1,j1,k1)=dpy+dpie+dpjs+dpkf
                    q03(i1,j1,k1+1)=dpy+dpiw+dpjs+dpkb
                    q03(i1+1,j1,k1+1)=dpy+dpie+dpjs+dpkb
                    q03(i1,j1+1,k1)=dpy+dpiw+dpjn+dpkf
                    q03(i1+1,j1+1,k1)=dpy+dpie+dpjn+dpkf
                    q03(i1,j1+1,k1+1)=dpy+dpiw+dpjn+dpkb
                    q03(i1+1,j1+1,k1+1)=dpy+dpie+dpjn+dpkb

                    dpy=py4(i,j,k)
                    dpiw=dp(py4(iw,j,k),dpy,py4(ie,j,k))
                    dpie=dp(py4(ie,j,k),dpy,py4(iw,j,k))
                    dpjs=dp(py4(i,js,k),dpy,py4(i,jn,k))
                    dpjn=dp(py4(i,jn,k),dpy,py4(i,js,k))
                    dpkf=dp(py4(i,j,kf),dpy,py4(i,j,kb))
                    dpkb=dp(py4(i,j,kb),dpy,py4(i,j,kf))
                    q04(i1,j1,k1)=dpy+dpiw+dpjs+dpkf
                    q04(i1+1,j1,k1)=dpy+dpie+dpjs+dpkf
                    q04(i1,j1,k1+1)=dpy+dpiw+dpjs+dpkb
                    q04(i1+1,j1,k1+1)=dpy+dpie+dpjs+dpkb
                    q04(i1,j1+1,k1)=dpy+dpiw+dpjn+dpkf
                    q04(i1+1,j1+1,k1)=dpy+dpie+dpjn+dpkf
                    q04(i1,j1+1,k1+1)=dpy+dpiw+dpjn+dpkb
                    q04(i1+1,j1+1,k1+1)=dpy+dpie+dpjn+dpkb

                    dpy=py5(i,j,k)
                    dpiw=dp(py5(iw,j,k),dpy,py5(ie,j,k))
                    dpie=dp(py5(ie,j,k),dpy,py5(iw,j,k))
                    dpjs=dp(py5(i,js,k),dpy,py5(i,jn,k))
                    dpjn=dp(py5(i,jn,k),dpy,py5(i,js,k))
                    dpkf=dp(py5(i,j,kf),dpy,py5(i,j,kb))
                    dpkb=dp(py5(i,j,kb),dpy,py5(i,j,kf))
                    q05(i1,j1,k1)=dpy+dpiw+dpjs+dpkf
                    q05(i1+1,j1,k1)=dpy+dpie+dpjs+dpkf
                    q05(i1,j1,k1+1)=dpy+dpiw+dpjs+dpkb
                    q05(i1+1,j1,k1+1)=dpy+dpie+dpjs+dpkb
                    q05(i1,j1+1,k1)=dpy+dpiw+dpjn+dpkf
                    q05(i1+1,j1+1,k1)=dpy+dpie+dpjn+dpkf
                    q05(i1,j1+1,k1+1)=dpy+dpiw+dpjn+dpkb
                    q05(i1+1,j1+1,k1+1)=dpy+dpie+dpjn+dpkb
                end do
            end do
        end do
        nx=nnx(ig+1)
        ny=nny(ig+1)
        nz=nnz(ig+1)
        py1(1:nx,1:ny,1:nz)=q01(1:nx,1:ny,1:nz)
        py2(1:nx,1:ny,1:nz)=q02(1:nx,1:ny,1:nz)
        py3(1:nx,1:ny,1:nz)=q03(1:nx,1:ny,1:nz)
        py4(1:nx,1:ny,1:nz)=q04(1:nx,1:ny,1:nz)
        py5(1:nx,1:ny,1:nz)=q05(1:nx,1:ny,1:nz)
    end do
    q31(1:nx,1:ny,1:nz,n)=q31(1:nx,1:ny,1:nz,n)+py1(1:nx,1:ny,1:nz)
    q32(1:nx,1:ny,1:nz,n)=q32(1:nx,1:ny,1:nz,n)+py2(1:nx,1:ny,1:nz)
    q33(1:nx,1:ny,1:nz,n)=q33(1:nx,1:ny,1:nz,n)+py3(1:nx,1:ny,1:nz)
    q34(1:nx,1:ny,1:nz,n)=q34(1:nx,1:ny,1:nz,n)+py4(1:nx,1:ny,1:nz)
    q35(1:nx,1:ny,1:nz,n)=q35(1:nx,1:ny,1:nz,n)+py5(1:nx,1:ny,1:nz)
    nx=nnx(ign)
    ny=nny(ign)
    nz=nnz(ign)
    end subroutine update

real function dp(a,b,c) !外部函数,处理当前点与左右点的关系
    implicit none
    real(8) :: a,b,c
    dp=(3.d0*a-2.d0*b-c)/64.d0
    return
    end

subroutine update2 !插值跳阶段
    use global
    implicit none
    integer ::ig,i1,j1,k1,iw,ie,jn,js,kf,kb
    real,external :: dp
    real(8) :: dpy,dpiw,dpie,dpjn,dpjs,dpkf,dpkb
    real(8),allocatable ::s1(:,:),s2(:,:),s3(:,:),s4(:,:),s5(:,:),s6(:,:)
    nx=nnx(nng)
    ny=nny(nng)
    nz=nnz(nng)
    do n=1,nt
        py1(1:nx,1:ny,1:nz)=q11(1:nx,1:ny,1:nz,n)
        py2(1:nx,1:ny,1:nz)=q12(1:nx,1:ny,1:nz,n)
        py3(1:nx,1:ny,1:nz)=q13(1:nx,1:ny,1:nz,n)
        py4(1:nx,1:ny,1:nz)=q14(1:nx,1:ny,1:nz,n)
        py5(1:nx,1:ny,1:nz)=q15(1:nx,1:ny,1:nz,n)
        py6(1:nx,1:ny,1:nz)=q16(1:nx,1:ny,1:nz,n)
        do k=1,nz
            k1=2*k-1
            kf=max(1,k-1)
            kb=min(nz,k+1)
            do j=1,ny
                j1=2*j-1
                js=max(1,j-1)
                jn=min(ny,j+1)
                do i=1,nx
                    i1=2*i-1
                    iw=max(1,i-1)
                    ie=min(nx,i+1)

                    dpy=py1(i,j,k)
                    dpiw=dp(py1(iw,j,k),dpy,py1(ie,j,k))
                    dpie=dp(py1(ie,j,k),dpy,py1(iw,j,k))
                    dpjs=dp(py1(i,js,k),dpy,py1(i,jn,k))
                    dpjn=dp(py1(i,jn,k),dpy,py1(i,js,k))
                    dpkf=dp(py1(i,j,kf),dpy,py1(i,j,kb))
                    dpkb=dp(py1(i,j,kb),dpy,py1(i,j,kf))
                    q11(i1,j1,k1,n)=dpy+dpiw+dpjs+dpkf
                    q11(i1+1,j1,k1,n)=dpy+dpie+dpjs+dpkf
                    q11(i1,j1,k1+1,n)=dpy+dpiw+dpjs+dpkb
                    q11(i1+1,j1,k1+1,n)=dpy+dpie+dpjs+dpkb
                    q11(i1,j1+1,k1,n)=dpy+dpiw+dpjn+dpkf
                    q11(i1+1,j1+1,k1,n)=dpy+dpie+dpjn+dpkf
                    q11(i1,j1+1,k1+1,n)=dpy+dpiw+dpjn+dpkb
                    q11(i1+1,j1+1,k1+1,n)=dpy+dpie+dpjn+dpkb

                    dpy=py2(i,j,k)
                    dpiw=dp(py2(iw,j,k),dpy,py2(ie,j,k))
                    dpie=dp(py2(ie,j,k),dpy,py2(iw,j,k))
                    dpjs=dp(py2(i,js,k),dpy,py2(i,jn,k))
                    dpjn=dp(py2(i,jn,k),dpy,py2(i,js,k))
                    dpkf=dp(py2(i,j,kf),dpy,py2(i,j,kb))
                    dpkb=dp(py2(i,j,kb),dpy,py2(i,j,kf))
                    q12(i1,j1,k1,n)=dpy+dpiw+dpjs+dpkf
                    q12(i1+1,j1,k1,n)=dpy+dpie+dpjs+dpkf
                    q12(i1,j1,k1+1,n)=dpy+dpiw+dpjs+dpkb
                    q12(i1+1,j1,k1+1,n)=dpy+dpie+dpjs+dpkb
                    q12(i1,j1+1,k1,n)=dpy+dpiw+dpjn+dpkf
                    q12(i1+1,j1+1,k1,n)=dpy+dpie+dpjn+dpkf
                    q12(i1,j1+1,k1+1,n)=dpy+dpiw+dpjn+dpkb
                    q12(i1+1,j1+1,k1+1,n)=dpy+dpie+dpjn+dpkb

                    dpy=py3(i,j,k)
                    dpiw=dp(py3(iw,j,k),dpy,py3(ie,j,k))
                    dpie=dp(py3(ie,j,k),dpy,py3(iw,j,k))
                    dpjs=dp(py3(i,js,k),dpy,py3(i,jn,k))
                    dpjn=dp(py3(i,jn,k),dpy,py3(i,js,k))
                    dpkf=dp(py3(i,j,kf),dpy,py3(i,j,kb))
                    dpkb=dp(py3(i,j,kb),dpy,py3(i,j,kf))
                    q13(i1,j1,k1,n)=dpy+dpiw+dpjs+dpkf
                    q13(i1+1,j1,k1,n)=dpy+dpie+dpjs+dpkf
                    q13(i1,j1,k1+1,n)=dpy+dpiw+dpjs+dpkb
                    q13(i1+1,j1,k1+1,n)=dpy+dpie+dpjs+dpkb
                    q13(i1,j1+1,k1,n)=dpy+dpiw+dpjn+dpkf
                    q13(i1+1,j1+1,k1,n)=dpy+dpie+dpjn+dpkf
                    q13(i1,j1+1,k1+1,n)=dpy+dpiw+dpjn+dpkb
                    q13(i1+1,j1+1,k1+1,n)=dpy+dpie+dpjn+dpkb

                    dpy=py4(i,j,k)
                    dpiw=dp(py4(iw,j,k),dpy,py4(ie,j,k))
                    dpie=dp(py4(ie,j,k),dpy,py4(iw,j,k))
                    dpjs=dp(py4(i,js,k),dpy,py4(i,jn,k))
                    dpjn=dp(py4(i,jn,k),dpy,py4(i,js,k))
                    dpkf=dp(py4(i,j,kf),dpy,py4(i,j,kb))
                    dpkb=dp(py4(i,j,kb),dpy,py4(i,j,kf))
                    q14(i1,j1,k1,n)=dpy+dpiw+dpjs+dpkf
                    q14(i1+1,j1,k1,n)=dpy+dpie+dpjs+dpkf
                    q14(i1,j1,k1+1,n)=dpy+dpiw+dpjs+dpkb
                    q14(i1+1,j1,k1+1,n)=dpy+dpie+dpjs+dpkb
                    q14(i1,j1+1,k1,n)=dpy+dpiw+dpjn+dpkf
                    q14(i1+1,j1+1,k1,n)=dpy+dpie+dpjn+dpkf
                    q14(i1,j1+1,k1+1,n)=dpy+dpiw+dpjn+dpkb
                    q14(i1+1,j1+1,k1+1,n)=dpy+dpie+dpjn+dpkb

                    dpy=py5(i,j,k)
                    dpiw=dp(py5(iw,j,k),dpy,py5(ie,j,k))
                    dpie=dp(py5(ie,j,k),dpy,py5(iw,j,k))
                    dpjs=dp(py5(i,js,k),dpy,py5(i,jn,k))
                    dpjn=dp(py5(i,jn,k),dpy,py5(i,js,k))
                    dpkf=dp(py5(i,j,kf),dpy,py5(i,j,kb))
                    dpkb=dp(py5(i,j,kb),dpy,py5(i,j,kf))
                    q15(i1,j1,k1,n)=dpy+dpiw+dpjs+dpkf
                    q15(i1+1,j1,k1,n)=dpy+dpie+dpjs+dpkf
                    q15(i1,j1,k1+1,n)=dpy+dpiw+dpjs+dpkb
                    q15(i1+1,j1,k1+1,n)=dpy+dpie+dpjs+dpkb
                    q15(i1,j1+1,k1,n)=dpy+dpiw+dpjn+dpkf
                    q15(i1+1,j1+1,k1,n)=dpy+dpie+dpjn+dpkf
                    q15(i1,j1+1,k1+1,n)=dpy+dpiw+dpjn+dpkb
                    q15(i1+1,j1+1,k1+1,n)=dpy+dpie+dpjn+dpkb

                    dpy=py6(i,j,k)
                    dpiw=dp(py6(iw,j,k),dpy,py6(ie,j,k))
                    dpie=dp(py6(ie,j,k),dpy,py6(iw,j,k))
                    dpjs=dp(py6(i,js,k),dpy,py6(i,jn,k))
                    dpjn=dp(py6(i,jn,k),dpy,py6(i,js,k))
                    dpkf=dp(py6(i,j,kf),dpy,py6(i,j,kb))
                    dpkb=dp(py6(i,j,kb),dpy,py6(i,j,kf))
                    q16(i1,j1,k1,n)=dpy+dpiw+dpjs+dpkf
                    q16(i1+1,j1,k1,n)=dpy+dpie+dpjs+dpkf
                    q16(i1,j1,k1+1,n)=dpy+dpiw+dpjs+dpkb
                    q16(i1+1,j1,k1+1,n)=dpy+dpie+dpjs+dpkb
                    q16(i1,j1+1,k1,n)=dpy+dpiw+dpjn+dpkf
                    q16(i1+1,j1+1,k1,n)=dpy+dpie+dpjn+dpkf
                    q16(i1,j1+1,k1+1,n)=dpy+dpiw+dpjn+dpkb
                    q16(i1+1,j1+1,k1+1,n)=dpy+dpie+dpjn+dpkb
                end do
            end do
        end do
    end do
    end subroutine update2

subroutine wwma  !求相对马赫数
    use global
    implicit none
    real(8) ::qqw,a

    do k=1,nz
        do j=1,ny
            do i=1,nx
                y1=yy0(i,j,k)
                z1=zz0(i,j,k)
                vx=pvx(i,j,k)
                vy=pvy(i,j,k)
                vz=pvz(i,j,k)
                wx=vx
                wy=vy+rpm*z1
                wz=vz-rpm*y1
                qqw=wx*wx+wy*wy+wz*wz
                a=1.4d0*p(i,j,k)/Q11(i,j,k,n)
                wma(i,j,k)=sqrt(qqw/a)
            end do
        end do
    end do
end subroutine wwma

 subroutine flow !进出口流量
    use global
    implicit none
    real(8) :: vf1,vf2,qin,qout,qinn,qoutt
    real(8),allocatable :: qin0(:),qout0(:)

            qin=0.d0
            qout=0.d0
            do n=1,nt
                do j=1,ny
                    do k=1,nz
                        vf1=-(s2x(1,j,k)*q12(1,j,k,n)+s2y(1,j,k)*q13(1,j,k,n)+s2z(1,j,k)*q14(1,j,k,n))
                        qin=qin+vf1

                        vf2=-(s2x(nx+1,j,k)*q12(nx,j,k,n)+s2y(nx+1,j,k)*q13(nx,j,k,n)+s2z(nx+1,j,k)*q14(nx,j,k,n))
                        qout=qout+vf2
                    end do
                end do
            end do
            qin=qin/dble(nt)
            qout=qout/dble(nt)
            if(myid==0)then
                allocate(qin0(0:lbb-1))
                allocate(qout0(0:lbb-1))
                qin0(0)=qin
                qinn=0.0
                qout0(0)=qout
                qoutt=0.0
                    do j=1,lbb-1
                        call MPI_RECV(qin0(j),1,MPI_DOUBLE_PRECISION,j,3,MPI_COMM_WORLD,status,ierr)
                        call MPI_RECV(qout0(j),1,MPI_DOUBLE_PRECISION,j,4,MPI_COMM_WORLD,status,ierr)
                    end do
                    do j=0,lbb-1
                        qinn=qinn+qin0(j)
                        qoutt=qoutt+qout0(j)
                    end do
                deallocate(qin0)
                deallocate(qout0)
                write(404+myid,"(i5,F15.6)")nitt,qinn
                write(707+myid,"(i5,F15.6)")nitt,qoutt
            else
                call MPI_SEND(qin,1,MPI_DOUBLE_PRECISION,0,3,MPI_COMM_WORLD,ierr)
                call MPI_SEND(qout,1,MPI_DOUBLE_PRECISION,0,4,MPI_COMM_WORLD,ierr)
            end if
    end subroutine flow

subroutine output  !输出
    use global
    implicit none
    integer :: kk,ifine,jfine,kfine,spa

        x(1:nx+1,1:ny+1,1:nz+1)=xf(1:nx+1,1:ny+1,1:nz+1)
        y(1:nx+1,1:ny+1,1:nz+1)=yf(1:nx+1,1:ny+1,1:nz+1)
        z(1:nx+1,1:ny+1,1:nz+1)=zf(1:nx+1,1:ny+1,1:nz+1)
        call in1out
        call lbout
        spa=30
        call span(spa)
        spa=95
        call span(spa)
    end subroutine output

subroutine in1out !进出口量,输出进口流量
    use global
    implicit none
    integer :: mml
    real(8) :: vf1,vf2,qin,qout1,qq2,tem,h0,t2,p2,pout,tout,d1,d2,eff,temp,ftt,fpp
    real(8),allocatable :: rr0(:),qout(:),fp(:),ft(:),qin0(:),qout0(:),d10(:),d20(:),eff0(:)
    real(8)::qinn,qoutt,d11,d22,efff

      if(myid==0)then
        open(4,file='inflow.dat') !时间谱配置点对精度的影响
        open(5,file='outflow.dat')
        open(8,file='pbi.dat')
        open(9,file='tbi.dat')
        open(11,file='eff.dat')
      end if
      do n=1,nt
          qin=0.d0
            do j=1,ny
                do k=1,nz
                    vf1=-(s2x(1,j,k)*q12(1,j,k,n)+s2y(1,j,k)*q13(1,j,k,n)+s2z(1,j,k)*q14(1,j,k,n))
                    qin=qin+vf1
                end do
            end do
            if(myid==0)then
                allocate(qin0(0:lbb-1))
                qin0(0)=qin
                qinn=0.0
                    do j=1,lbb-1
                        call MPI_RECV(qin0(j),1,MPI_DOUBLE_PRECISION,j,3,MPI_COMM_WORLD,status,ierr)
                    end do
                    do j=0,lbb-1
                        qinn=qinn+qin0(j)
                    end do
                deallocate(qin0)
                write(4,"(f6.3,F15.6)")dble(n-1)/dble(nt),qinn
            else
                call MPI_SEND(qin,1,MPI_DOUBLE_PRECISION,0,3,MPI_COMM_WORLD,ierr)
            end if
      end do

    allocate(rr0(ny))
    allocate(qout(ny))
    allocate(fp(ny))
    allocate(ft(ny))
      do j=1,ny
          y1=yy0(nx,j,1)
          z1=zz0(nx,j,1)
          rr0(j)=sqrt(y1*y1+z1*z1)
      end do
      do n=1,nt
          qout=0.d0
          fp=0.d0
          ft=0.d0
          qout1=0.d0
          fpp=0.d0
          ftt=0.d0
          write(zonename,"(i3)") n
          open(44,file='pbispan-timel'//trim(adjustl(zonename))//'.dat')
          open(55,file='tbispan-timel'//trim(adjustl(zonename))//'.dat')
          open(66,file='effspan-timel'//trim(adjustl(zonename))//'.dat')
          do j=1,ny
             temp=(rr0(j)-rr0(1))/(rr0(ny)-rr0(1))*100.d0
             do k=1,nz
                 vf2=-(s2x(nx+1,j,k)*q12(nx,j,k,n)+s2y(nx+1,j,k)*q13(nx,j,k,n)+s2z(nx+1,j,k)*q14(nx,j,k,n))
                 vx=q12(nx,j,k,n)/q11(nx,j,k,n)
                 vy=q13(nx,j,k,n)/q11(nx,j,k,n)
                 vz=q14(nx,j,k,n)/q11(nx,j,k,n)
                 qq2=vx*vx+vy*vy+vz*vz
                 pp=0.4d0*(q15(nx,j,k,n)-0.5d0*q11(nx,j,k,n)*qq2)
                 h0=(q15(nx,j,k,n)+pp)/q11(nx,j,k,n)
                 t2=h0/cp
                 p2=pp*(1.d0-0.5d0*qq2/h0)**(-3.5)
                 ft(j)=ft(j)+t2*vf2
                 fp(j)=fp(j)+p2*vf2
                 qout(j)=qout(j)+vf2
             end do
             pout=fp(j)/qout(j)
             tout=ft(j)/qout(j)
             d1=pout/pt
             d2=tout/(ht/cp)
             eff=(d1**(2.d0/7.d0)-1.d0)/(d2-1.d0)*100.d0
            if(myid==0)then
                allocate(d10(0:lbb-1))
                allocate(d20(0:lbb-1))
                allocate(eff0(0:lbb-1))
                d10(0)=d1
                d20(0)=d2
                eff0(0)=eff
                d11=0.0
                d22=0.0
                efff=0.0
                    do i=1,lbb-1
                        call MPI_RECV(d10(i),1,MPI_DOUBLE_PRECISION,i,3,MPI_COMM_WORLD,status,ierr)
                        call MPI_RECV(d20(i),1,MPI_DOUBLE_PRECISION,i,4,MPI_COMM_WORLD,status,ierr)
                        call MPI_RECV(eff0(i),1,MPI_DOUBLE_PRECISION,i,5,MPI_COMM_WORLD,status,ierr)
                    end do
                    do i=0,lbb-1
                        d11=d11+d10(i)
                        d22=d22+d20(i)
                        efff=efff+eff0(i)
                    end do
                deallocate(d10)
                deallocate(d20)
                deallocate(eff0)
                write(44,"(2F15.6)")d11/(mml*lbb),temp
                write(55,"(2F15.6)")d22/(mml*lbb),temp
                write(66,"(2F15.6)")efff/(mml*lbb),temp
            else
                call MPI_SEND(d1,1,MPI_DOUBLE_PRECISION,0,3,MPI_COMM_WORLD,ierr)
                call MPI_SEND(d2,1,MPI_DOUBLE_PRECISION,0,4,MPI_COMM_WORLD,ierr)
                call MPI_SEND(eff,1,MPI_DOUBLE_PRECISION,0,5,MPI_COMM_WORLD,ierr)
            end if
             ftt=ftt+ft(j)
             fpp=fpp+fp(j)
             qout1=qout1+qout(j)
          end do
          close(44)
          close(55)
          close(66)
          pout=fpp/qout1
          tout=ftt/qout1
          d1=pout/pt
          d2=tout/(ht/cp)
          eff=(d1**(2.d0/7.d0)-1.d0)/(d2-1.d0)*100.d0
            if(myid==0)then
                allocate(qout0(0:lbb-1))
                allocate(d10(0:lbb-1))
                allocate(d20(0:lbb-1))
                allocate(eff0(0:lbb-1))
                qout0(0)=qout1
                d10(0)=d1
                d20(0)=d2
                eff0(0)=eff
                qoutt=0.0
                d11=0.0
                d22=0.0
                efff=0.0
                    do j=1,lbb-1
                        call MPI_RECV(qout0(j),1,MPI_DOUBLE_PRECISION,j,2,MPI_COMM_WORLD,status,ierr)
                        call MPI_RECV(d10(j),1,MPI_DOUBLE_PRECISION,j,3,MPI_COMM_WORLD,status,ierr)
                        call MPI_RECV(d20(j),1,MPI_DOUBLE_PRECISION,j,4,MPI_COMM_WORLD,status,ierr)
                        call MPI_RECV(eff0(j),1,MPI_DOUBLE_PRECISION,j,5,MPI_COMM_WORLD,status,ierr)
                    end do
                    do j=0,lbb-1
                        qoutt=qoutt+qout0(j)
                        d11=d11+d10(j)
                        d22=d22+d20(j)
                        efff=efff+eff0(j)
                    end do
                deallocate(qout0)
                deallocate(d10)
                deallocate(d20)
                deallocate(eff0)
                write(5,"(f6.3,F15.6)")dble(n-1)/dble(nt),qoutt
                write(8,"(f6.3,F15.6)")dble(n-1)/dble(nt),d11/lbb
                write(9,"(f6.3,F15.6)")dble(n-1)/dble(nt),d22/lbb
                write(11,"(f6.3,F15.6,A)")dble(n-1)/dble(nt),efff/lbb
            else
                call MPI_SEND(qout1,1,MPI_DOUBLE_PRECISION,0,2,MPI_COMM_WORLD,ierr)
                call MPI_SEND(d1,1,MPI_DOUBLE_PRECISION,0,3,MPI_COMM_WORLD,ierr)
                call MPI_SEND(d2,1,MPI_DOUBLE_PRECISION,0,4,MPI_COMM_WORLD,ierr)
                call MPI_SEND(eff,1,MPI_DOUBLE_PRECISION,0,5,MPI_COMM_WORLD,ierr)
            end if
      end do
    deallocate(qout)
    deallocate(fp)
    deallocate(ft)
    deallocate(rr0)
    end subroutine in1out

subroutine lbout   !叶片 输出叶片计算结果
    use global
    implicit none
    integer :: e,lmm,h
    real(8) ::temp,t1,t2,cor1,sir1,qq2,cvl,qqw,a

    open(44,file='out-'//trim(adjustl(id_m))//'myid.dat')
    write(44,*) "VARIABLES=x,y,z,u,v,w,wma,pressure,density,ut"
    do n=1,nt
        call ppp
        call wwma
        write(zonename,"(i3)") n
        zonename='out_'//trim(adjustl(zonename))
        if(n==1) then
            write(44,"(1x,2A,2x,A,i4,2x,A,i4,2x,A,i4,2x,A,F8.4)") "ZONE T=",zonename, "I=",nx+1,"J=",ny+1,"K=",nz+1,"DATAPACKING=BLOCK, VARLOCATION=([4-10]=CELLCENTERED),SOLUTIONTIME=",dble(n-1)/dble(nt)
            write(44,*) x(1:nx+1,1:ny+1,1:nz+1)
            write(44,*) y(1:nx+1,1:ny+1,1:nz+1)
            write(44,*) z(1:nx+1,1:ny+1,1:nz+1)
            write(44,*) pvx(1:nx,1:ny,1:nz)
            write(44,*) pvy(1:nx,1:ny,1:nz)
            write(44,*) pvz(1:nx,1:ny,1:nz)
            write(44,*) wma(1:nx,1:ny,1:nz)
            write(44,*) p(1:nx,1:ny,1:nz)
            write(44,*) q11(1:nx,1:ny,1:nz,n)
            write(44,*) q16(1:nx,1:ny,1:nz,n)
        else
            write(44,"(1x,2A,2x,A,i4,2x,A,i4,2x,A,i4,2x,A,F8.4)") "ZONE T=",zonename, "I=",nx+1,"J=",ny+1,"K=",nz+1,"DATAPACKING=BLOCK,VARSHARELIST=([1-3]=1), VARLOCATION=([4-10]=CELLCENTERED),SOLUTIONTIME=",dble(n-1)/dble(nt)
            write(44,*) pvx(1:nx,1:ny,1:nz)
            write(44,*) pvy(1:nx,1:ny,1:nz)
            write(44,*) pvz(1:nx,1:ny,1:nz)
            write(44,*) wma(1:nx,1:ny,1:nz)
            write(44,*) p(1:nx,1:ny,1:nz)
            write(44,*) q11(1:nx,1:ny,1:nz,n)
            write(44,*) q16(1:nx,1:ny,1:nz,n)
        end if
    end do
    close(44)
    end subroutine lbout

subroutine span(spa) !叶高 输出沿叶高结果
    use global
    implicit none
    integer ::e
    integer,intent(in) ::spa
    integer ::h,lblb,xll,zll
    real(8) :: temp,tem,cvl,x1,x2,x3,r1,t1,t2,qq2,qqw,qqr,a,ht0,ppt,cor1,sir1
    real(8),allocatable :: hx(:),hy(:),hz(:),hr(:),s(:),hu(:),hv(:),hw(:),hp(:),hpt(:),hh(:),hmiu(:),hwma(:)
    real(8),allocatable :: hxx(:,:),hyy(:,:),hzz(:,:),huu(:,:,:),hvv(:,:,:),hww(:,:,:),hpp(:,:,:),hppt(:,:,:),hht(:,:,:),hwwma(:,:,:),hmmiu(:,:,:)

    allocate(hx(1:ny+1))
    allocate(hy(1:ny+1))
    allocate(hz(1:ny+1))
    allocate(hr(1:ny+1))
    allocate(s(1:ny+1))
    allocate(hu(1:ny))
    allocate(hv(1:ny))
    allocate(hw(1:ny))
    allocate(hp(1:ny))
    allocate(hpt(1:ny))
    allocate(hh(1:ny))
    allocate(hmiu(1:ny))
    allocate(hwma(1:ny))
    allocate(hxx(1:nx+1,1:nz+1))
    allocate(hyy(1:nx+1,1:nz+1))
    allocate(hzz(1:nx+1,1:nz+1))
    allocate(huu(1:nx,1:nz,1:nt))
    allocate(hvv(1:nx,1:nz,1:nt))
    allocate(hww(1:nx,1:nz,1:nt))
    allocate(hpp(1:nx,1:nz,1:nt))
    allocate(hppt(1:nx,1:nz,1:nt))
    allocate(hht(1:nx,1:nz,1:nt))
    allocate(hwwma(1:nx,1:nz,1:nt))
    allocate(hmmiu(1:nx,1:nz,1:nt))

    temp=dble(spa)/100.d0
    do k=1,nz+1
        do i=1,nx+1
            s(1)=0.d0
            do j=1,ny+1
                hx(j)=x(i,j,k)
                hy(j)=y(i,j,k)
                hz(j)=z(i,j,k)
                hr(j)=sqrt(hy(j)*hy(j)+hz(j)*hz(j))
                if(j>1)then
                    s(j)=s(j-1)+sqrt((hx(j)-hx(j-1))**2+(hr(j)-hr(j-1))**2)
                end if
            end do
            t1=s(ny+1)*temp
            call wl(s(1:ny+1),hx(1:ny+1),t1,hxx(i,k),ny+1,1)
            call wl(s(1:ny+1),hy(1:ny+1),t1,hyy(i,k),ny+1,1)
            call wl(s(1:ny+1),hz(1:ny+1),t1,hzz(i,k),ny+1,1)
        end do
    end do
    do k=1,nz
        do i=1,nx
            s(1)=0.d0
            do j=1,ny
                hx(j)=xx0(i,j,k)
                hy(j)=yy0(i,j,k)
                hz(j)=zz0(i,j,k)
                hr(j)=sqrt(hy(j)*hy(j)+hz(j)*hz(j))
                if(j>1)then
                    s(j)=s(j-1)+sqrt((hx(j)-hx(j-1))**2+(hr(j)-hr(j-1))**2)
                end if
            end do
            t1=s(ny)*temp
            do n=1,nt
                do j=1,ny
                    y1=yy0(i,j,k)
                    z1=zz0(i,j,k)
                    hu(j)=Q12(i,j,k,n)/Q11(i,j,k,n)
                    hv(j)=Q13(i,j,k,n)/Q11(i,j,k,n)
                    hw(j)=Q14(i,j,k,n)/Q11(i,j,k,n)
                    wx=hu(j)
                    wy=hv(j)+rpm*z1
                    wz=hw(j)-rpm*y1
                    qqw=wx*wx+wy*wy+wz*wz
                    qqr=rpm*rpm*(z1*z1+y1*y1)
                    qq2=hu(j)*hu(j)+hv(j)*hv(j)+hw(j)*hw(j)
                    hp(j)=0.4d0*(Q15(i,j,k,n)-0.5d0*Q11(i,j,k,n)*qq2)
                    tem=hp(j)/(Q11(i,j,k,n)*rg)
                    cvl=cvl0*((tem/t0)**1.5)*(t0+ts)/(tem+ts)
                    a=1.4d0*hp(j)/Q11(i,j,k,n)
                    hwma(j)=sqrt(qqw/a)
                    hh(j)=(Q15(i,j,k,n)+hp(j))/Q11(i,j,k,n)
                    hpt(j)=hp(j)*(1.d0-0.5d0*qq2/hh(j))**(-3.5)
                    hmiu(j)=Q16(i,j,k,n)/cvl
                end do
                call wl(s(1:ny),hu(1:ny),t1,huu(i,k,n),ny,1)
                call wl(s(1:ny),hv(1:ny),t1,hvv(i,k,n),ny,1)
                call wl(s(1:ny),hw(1:ny),t1,hww(i,k,n),ny,1)
                call wl(s(1:ny),hp(1:ny),t1,hpp(i,k,n),ny,1)
                call wl(s(1:ny),hpt(1:ny),t1,hppt(i,k,n),ny,1)
                call wl(s(1:ny),hh(1:ny),t1,hht(i,k,n),ny,1)
                call wl(s(1:ny),hwma(1:ny),t1,hwwma(i,k,n),ny,1)
                call wl(s(1:ny),hmiu(1:ny),t1,hmmiu(i,k,n),ny,1)
            end do
        end do
    end do
    !*****************输出
    write(nnspan,'(i3)') spa
    open(21,file=trim(adjustl(nnspan))//'%span-'//trim(adjustl(id_m))//'myid.dat')
    write(21,*) "VARIABLES=x,y,z,u,v,w,wma,pressure,pt,zht,ut"
    do n=1,nt
        write(zonename,"(i3)") n
        zonename='out_'//trim(adjustl(zonename))
        if(n==1) then
            write(21,"(1x,2A,2x,A,i4,2x,A,i4,2x,A,i4,2x,A,F8.4)") "ZONE T=",zonename,"I=",nx+1,"J=",1,"K=",nz+1,"DATAPACKING=BLOCK, VARLOCATION=([4-11]=CELLCENTERED),SOLUTIONTIME=",dble(n-1)/dble(nt)
            write(21,*) hxx(1:nx+1,1:nz+1)
            write(21,*) hyy(1:nx+1,1:nz+1)
            write(21,*) hzz(1:nx+1,1:nz+1)
            write(21,*) huu(1:nx,1:nz,n)
            write(21,*) hvv(1:nx,1:nz,n)
            write(21,*) hww(1:nx,1:nz,n)
            write(21,*) hwwma(1:nx,1:nz,n)
            write(21,*) hpp(1:nx,1:nz,n)
            write(21,*) hppt(1:nx,1:nz,n)
            write(21,*) hht(1:nx,1:nz,n)
            write(21,*) hmmiu(1:nx,1:nz,n)
        else
            write(21,"(1x,2A,2x,A,i4,2x,A,i4,2x,A,i4,2x,A,F8.4)") "ZONE T=",zonename,"I=",nx+1,"J=",1,"K=",nz+1,"DATAPACKING=BLOCK,VARSHARELIST=([1-3]=1), VARLOCATION=([4-11]=CELLCENTERED),SOLUTIONTIME=",dble(n-1)/dble(nt)
            write(21,*) huu(1:nx,1:nz,n)
            write(21,*) hvv(1:nx,1:nz,n)
            write(21,*) hww(1:nx,1:nz,n)
            write(21,*) hwwma(1:nx,1:nz,n)
            write(21,*) hpp(1:nx,1:nz,n)
            write(21,*) hppt(1:nx,1:nz,n)
            write(21,*) hht(1:nx,1:nz,n)
            write(21,*) hmmiu(1:nx,1:nz,n)
        end if
    end do
    close(21)
    end subroutine span

subroutine wl(medx,medr,spax,spar,n1,n2)    !三次样条插值
    implicit none
    integer :: i,j
    integer,intent(in) ::n1,n2
    real(8),intent(in) :: medx(n1),medr(n1),spax(n2)
    real(8),intent(out) :: spar(n2)
    real(8) ::h(1:n1-1),g(2:n1-1),lamt(2:n1-1),miu(2:n1-1),a(2:n1-1),m(1:n1)
    real(8) :: h1

    spar=0.d0
    m(1)=(medr(2)-medr(1))/(medx(2)-medx(1))
    m(n1)=(medr(n1)-medr(n1-1))/(medx(n1)-medx(n1-1))
    do i=1,n1-1
        h(i)=medx(i+1)-medx(i)
    end do
    do i=2,n1-1
        a(i)=2.d0
        lamt(i)=h(i)/(h(i)+h(i-1))
        miu(i)=1d0-lamt(i)
        g(i)=3.d0*(lamt(i)*(medr(i)-medr(i-1))/h(i-1)+miu(i)*(medr(i+1)-medr(i))/h(i))
    end do
    g(2)=g(2)-lamt(2)*m(1)
    g(n1-1)=g(n1-1)-miu(n1-1)*m(n1)
    call tdma(lamt(3:n1-1),a(2:n1-1),miu(2:n1-2),g(2:n1-1),m(2:n1-1),n1-2)
    do j=1,n2
        if(spax(j)<medx(n1))then
            i=1
            do while(spax(j)>medx(i+1))
                i=i+1
            end do
        else
            i=n1-1
        end if
        h1=(spax(j)-medx(i))/h(i)
        spar(j)=spar(j)+medr(i)*(2.d0*h1+1.d0)*(h1-1.d0)*(h1-1.d0)+m(i)*h(i)*h1*(h1-1.d0)*(h1-1.d0)
        h1=(medx(i+1)-spax(j))/h(i)
        spar(j)=spar(j)+medr(i+1)*(2.d0*h1+1.d0)*(h1-1.d0)*(h1-1.d0)-m(i+1)*h(i)*h1*(h1-1.d0)*(h1-1.d0)
    end do
end subroutine wl
