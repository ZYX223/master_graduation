Module global
    use mpi
    implicit none
    save
    integer ::m,n,mm,nn,mt,nt,nitt,ci,nci,nb,i,j,k,l,ii,jj,kk,ifine,jfine,kfine,nx,ny,nz,nx1,nx2,nxyp,ib,it,jb,jt,jt0,ny0,iread,cc,xl,yl,zl,xll,xlln,xlnn,yln,zln!最密网格点数，时间层数，多重网格数
    integer,allocatable ::lb(:),xln(:)
    real(8) ::y1,z1,rr,sir,cor,vx,vy,vz,wx,wy,wz,dim,en,pp,delt,cvl,a,qq2,two3,spa1,spa2
    real(8),parameter ::pi=3.141592654
    real(8) ::c2,cfl,a4,a2,beta1,beta2,aa1,aa2,aa3,val1,val2,val3,val4,val5,val6,val7,val8,val,dimin,ppin
    character(len=50):: nnspan,id_mm,id_sa
    real(8),allocatable ::Q01(:,:,:),Q02(:,:,:),Q03(:,:,:),Q04(:,:,:),Q05(:,:,:),Q06(:,:,:),gdf(:,:,:)   !RK循环前的守恒变量
    real(8),allocatable ::Q11(:,:,:),Q12(:,:,:),Q13(:,:,:),Q14(:,:,:),Q15(:,:,:),Q16(:,:,:) !守恒变量
    real(8),allocatable ::Q21(:,:,:,:),Q22(:,:,:,:),Q23(:,:,:,:),Q24(:,:,:,:),Q25(:,:,:,:),Q26(:,:,:,:) !保存前两时层量，二阶向后差分用到
    real(8),allocatable ::ts1(:,:,:),ts2(:,:,:),ts3(:,:,:),ts4(:,:,:),ts5(:,:,:),ts6(:,:,:) !二阶向后差分项的后两项
    real(8),allocatable ::AV1(:,:,:),AV2(:,:,:),AV3(:,:,:),AV4(:,:,:),AV5(:,:,:),AV6(:,:,:) !人工粘性项
    real(8),allocatable ::qc1(:,:,:),qc2(:,:,:),qc3(:,:,:),qc4(:,:,:),qc5(:,:,:),qc6(:,:,:) !对流通量
    real(8),allocatable ::qv1(:,:,:),qv2(:,:,:),qv3(:,:,:),qv4(:,:,:),qv5(:,:,:),qv6(:,:,:) !粘性通量
    real(8),allocatable ::x(:,:,:),y(:,:,:),z(:,:,:),x0(:,:,:),y0(:,:,:),z0(:,:,:),betax(:,:),betay(:,:),betaz(:,:),peb(:),xf0(:,:,:),yf0(:,:,:),zf0(:,:,:)
    real(8),allocatable ::xx01(:,:,:),yy01(:,:,:),zz01(:,:,:),xx02(:,:,:),yy02(:,:,:),zz02(:,:,:),xx03(:,:,:),yy03(:,:,:),zz03(:,:,:),xx0(:,:,:),yy0(:,:,:),zz0(:,:,:)
    real(8),allocatable ::s1x(:,:,:),s1y(:,:,:),s1z(:,:,:),s2x(:,:,:),s2y(:,:,:),s2z(:,:,:),s3x(:,:,:),s3y(:,:,:),s3z(:,:,:),vv0(:,:,:)
    real(8),allocatable ::pvx(:,:,:),pvy(:,:,:),pvz(:,:,:),vth(:,:,:),vre(:,:,:),p(:,:,:),t(:,:,:),time(:,:,:),temp0(:)
    real(8),allocatable ::rpm(:),ma(:),dmini(:,:,:),sri(:,:,:),srj(:,:,:),srk(:,:,:)
    real(8),allocatable :: hxx(:,:),hyy(:,:),hzz(:,:),spa0(:,:),hxx1(:,:),hyy1(:,:),hzz1(:,:),hxx2(:,:),hyy2(:,:),hzz2(:,:),spa01(:,:),spa02(:,:),hspa(:,:)
    real(8) ::ta,timl,pt,ht,rout,pb0,pb1,period,rmsm !总压，总焓，周期，动叶前后缘、叶根顶位置
    real(8) ::cvl0,t0,ts,cp,ccv,prt,prl,rg,cv1,cv2,kap,sigmav,cb1,cb2,cw1,cw2,cw3,cr1,cr2,cr3
    !并行量
    character(len=9):: id_m,id_rm,id_lm,secname,id_num1!进程号,叶排号,各叶排内叶片号
    real(8),allocatable :: xxs(:,:,:,:),yys(:,:,:,:),zzs(:,:,:,:),vv(:,:),v(:,:,:,:,:,:)
    integer :: cli(3),bk(3),ek(3),bkk(3),ekk(3)
    integer :: lm,lbm,slidm,rm,ssum,myid,myidl,myidm,myidr,numprocs,numpp,ierr,rc,status(MPI_STATUS_SIZE)
End module

program main
    use global
    implicit none
     
    call MPI_INIT( ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
    call ini           !主进程读入控制参数、全局网格坐标
    call allocation    !分配内存
    call ini2          !给1时层初始化
    call distance      !最细网格层上，求各网格中心点到固壁面最短的距离
    call geo           !各几何参数
    call spaa           !由于y方向分区，为后处理输出打基础
    call cfd         !spa1和spa2叶高截面物理量
    call MPI_FINALIZE(rc)
end
    
subroutine cfd!采用双时间步推进方法，按周期循环-物理时间循环-虚拟时间迭代循环来
    use global
    implicit none
    
     Periodcycle: do m=mm,mt       !周期数外循环
         if(m>mm)then
            nn=1
         end if
         Physicaltimelevel: do n=nn,nt !用到的物理时步中循环
             if(slidm/=0)then!动静叶网格线并不平行，故不能简单的面积比代表体积比，用三维！！
                 call slid          !动静叶传递，各通道得到相邻叶排所有通道的边界网格顶点坐标值
                 call overlapzone   !求出真实进程被哪些相邻进程裁剪，及裁剪的k在哪部分
                 call overlapgrid   !当前进程每个网格被拆分的面积比
                 call overlapadj    !重叠部分占当前进程范围，或者说《相邻进程>的重叠范围，后续bc发收物理量的范围，并传递交界面第一层网格x间距
             end if
             call tsd           !二阶向后差分的前两时层值
             Pseudotimelevel: do ci=1,nci    !用到的虚拟时步内循环
                 nitt=nitt+1
                 call march    !物理量计算1次又1步(4步)
                 call residual  !在最密网格上计算密度残差，并作为最后收敛的判断条件
                 write(20,"(i10,2x,F10.5,x,i3,x,i3,x,i5)") nitt,rmsm,ci,n,m !输出每一步迭代后的各时间层的最大残差，并计录每步的墙上时间
                 if(rmsm.gt.3d0) then
                 else if(rmsm.le.3d0) then
                 else
                     stop!*****判断残差若为NAN，则停止运行
                 endif
             end do Pseudotimelevel
             call store        !存储当前及前一步虚拟迭代值，并保存中间值
             call probe
             if(slidm/=0)then
                 call test
             end if
         end do Physicaltimelevel
         call output0
     end do Periodcycle
end subroutine cfd
    
subroutine ini  !主进程读入控制参数和坐标
    use global
    implicit none
    real(8) ::temp,t1,t2
    integer :: nxm1,nym1,nzm1,num1,ml
    integer,allocatable :: nxm(:),nym(:),nzm(:),ibm(:),itm(:),jbm(:),jtm(:)
    !***湍流计算所需参数
    cvl0=1.7161D-5
    t0=273.16d0!288.15
    ts=110.4d0 ! 124
    rg=287.d0
    cp=rg*1.4d0/0.4d0
    ccv=cp-rg
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
    cr2=2.d0!09文献为2，00文献为12
    cr3=1.d0
    two3=2./3.
     open(10,file="ini3ji.dat")       !动叶入口总焓，总压，中间级压力pb0，出口静压pb1
     read(10,*) iread,nt,beta1,beta2,cfl,a2,a4,ht,pt,pb1,c2,nb,mt,nci,spa1,spa2!c2,入口涡粘度因子。mt所用周期个数,nci是某一时刻下虚拟时步个数
     allocate(nxm(nb))
     allocate(nym(nb))
     allocate(nzm(nb))
     allocate(ibm(nb))
     allocate(itm(nb))
     allocate(jbm(nb))
     allocate(jtm(nb))
     allocate(lb(nb))
     allocate(rpm(nb))
     allocate(ma(nb))
     do l=1,nb
        read(10,*)
        read(10,*) nxm(l),nym(l),nzm(l),ibm(l),itm(l),jbm(l),jtm(l),lb(l),rpm(l),ma(l)   !itm,jtm代表固壁后的第一个流场网格，所以要减1
        rpm(l)=rpm(l)*pi/30.d0     !单位换算成rad/s--旋转角速度，决定周期
     end do
     close(10)
     period=2.d0*pi/abs(rpm(1))/lb(1)   !除动叶叶片数即为单通道（目前感觉除哪个叶排的叶片数都差不多）
     delt=period/dble(nt)!真实物理时间步长
     
     numpp=10
     if(myid<numpp*lb(1))then
         rm=1    !各进程所在叶排号1,2
     else
         rm=2
     end if
     write(id_m,'(i5)')myid
     write(id_rm,'(i3)') rm
     nx=nxm(rm)
     ny=nym(rm)
     nz=nzm(rm)
     ib=ibm(rm)
     it=itm(rm)-1
     jb=jbm(rm)
     jt=jtm(rm)-1
     jt0=jt
     allocate(xf0(1:nx+1,1:ny+1,1:nz+1))
     allocate(yf0(1:nx+1,1:ny+1,1:nz+1))
     allocate(zf0(1:nx+1,1:ny+1,1:nz+1))
     open(55,file='grid-'//trim(adjustl(id_rm))//'.dat')
     read(55,*) nxm1,nym1,nzm1
     do k=1,nzm1
         do j=1,nym1
             do i=1,nxm1
                 read(55,*)xf0(i,j,k),yf0(i,j,k),zf0(i,j,k)  !主进程读入总区坐标
             end do
         end do
     end do
     close(55)
     !各叶道分块
     xlnn=3      !x方向分三个大块的总数
     allocate(xln(0:xlnn-1))
     if(rm==1)then
         xln(0)=2!x方向每个大块的内部总数
         xln(2)=1
     else
         xln(0)=1
         xln(2)=2
     end if
     xln(1)=2
     yln=2
     zln=1
     num1=mod(myid,numpp)   !每个叶道内的进程编号0~79
     if(num1<xln(0)*yln*zln)then
         xlln=num1/(xln(0)*yln*zln)!x方向三大块的粗略编号0，1，2
     else
         xlln=(num1-xln(0)*yln*zln)/(xln(1)*yln*zln)+1
     end if
     ml=mod(num1,xln(xlln)*zln)  !在前两块内每个y上的8个进程，0~7，第三块有4个进程，0~3
     xl=ml/zln   !各大块的内部编号0~3，0~3，0~1
     xll=xl!综合三个大块下，各份在x方向的全局编号：0~9
     do i=1,xlln
         xll=xll+xln(i-1)
     end do
     do i=1,xlln
         num1=num1-xln(i-1)*yln*zln
     end do
     num1=mod(num1,xln(xlln)*yln*zln)!每个大块内的进程编号0~31，0~31，0~15
     yl=num1/(xln(xlln)*zln)  !y方向编号
     zl=mod(myid,zln)  !z方向编号
     !各分区坐标值
     nx=(it-ib+1)/xln(1) !当前进程的网格数
     if(rm==1)then
         nx1=nx!当前进程前的每进程网格数
         if(xll==xln(0)+xln(1)+xln(2)-1)then
             nx=24
         end if
     else
         nx1=24!当前进程前的每进程网格数
         if(xll==0)then
             nx=24
         end if
     end if
     ny=ny/yln
     nz=nz/zln
     allocate(x0(0:nx+2,1:ny+1,1:nz+1))
     allocate(y0(0:nx+2,1:ny+1,1:nz+1))
     allocate(z0(0:nx+2,1:ny+1,1:nz+1))
     if(rm==1)then
         x0(1:nx+1,1:ny+1,1:nz+1)=xf0(xll*nx1+1:xll*nx1+nx+1,yl*ny+1:(yl+1)*ny+1,zl*nz+1:(zl+1)*nz+1)
         y0(1:nx+1,1:ny+1,1:nz+1)=yf0(xll*nx1+1:xll*nx1+nx+1,yl*ny+1:(yl+1)*ny+1,zl*nz+1:(zl+1)*nz+1)
         z0(1:nx+1,1:ny+1,1:nz+1)=zf0(xll*nx1+1:xll*nx1+nx+1,yl*ny+1:(yl+1)*ny+1,zl*nz+1:(zl+1)*nz+1)
     else
         x0(1:nx+1,1:ny+1,1:nz+1)=xf0((xll-1)*nx+nx1+1:xll*nx+nx1+1,yl*ny+1:(yl+1)*ny+1,zl*nz+1:(zl+1)*nz+1)
         y0(1:nx+1,1:ny+1,1:nz+1)=yf0((xll-1)*nx+nx1+1:xll*nx+nx1+1,yl*ny+1:(yl+1)*ny+1,zl*nz+1:(zl+1)*nz+1)
         z0(1:nx+1,1:ny+1,1:nz+1)=zf0((xll-1)*nx+nx1+1:xll*nx+nx1+1,yl*ny+1:(yl+1)*ny+1,zl*nz+1:(zl+1)*nz+1)
     end if
     !旋转复制得整圈网格
     allocate(x(0:nx+2,1:ny+1,1:nz+1))
     allocate(y(0:nx+2,1:ny+1,1:nz+1))
     allocate(z(0:nx+2,1:ny+1,1:nz+1))
     lbm=myid/numpp!各进程所在的叶片号0,1,,35,36,,81
     if(lbm<lb(1))then
         lm=lbm
     else
         lm=lbm-lb(1)!经处理后，各进程所在叶排的叶片号，0,,35.0,,45
     end if
     temp=2.d0/dble(lb(rm))*pi*lm
     cor=cos(temp)                                       !0             ^
     sir=sin(temp)                                       !1             ^        叶片往上旋转
     do k=1,nz+1                                         !nz            ^
        do j=1,ny+1                                      !nz+1          1
            do i=1,nx+1
                x(i,j,k)=x0(i,j,k)
                y(i,j,k)=y0(i,j,k)*cor+z0(i,j,k)*sir !按叶片旋转方向旋转,若相反则是为了求0、nz+1坐标，1-nz的顺序与叶片旋转方向相反
                z(i,j,k)=-y0(i,j,k)*sir+z0(i,j,k)*cor
            end do
        end do
     end do
     !叶片在y方向范围
     if(yln==2)then
         if(yl==0)then
             jt=ny
         else
             jt=jt-ny
         end if
     end if
     !判断是否执行滑移网格
     slidm=0
     if(rm==1)then
         if(xll==xln(0)+xln(1)+xln(2)-1)then
             slidm=1
         end if
     else
         if(xll==0)then
             slidm=2
         end if
     end if
     deallocate(nxm)
     deallocate(nym)
     deallocate(nzm)
     deallocate(ibm)
     deallocate(itm)
     deallocate(jbm)
     deallocate(jtm)
    end subroutine ini

subroutine allocation  !申请内存
    use global
    implicit none

    allocate(gdf(15,0:ny+1,0:nz+1))
    allocate(q11(0:nx+1,0:ny+1,0:nz+1))
    allocate(q12(0:nx+1,0:ny+1,0:nz+1))
    allocate(q13(0:nx+1,0:ny+1,0:nz+1))
    allocate(q14(0:nx+1,0:ny+1,0:nz+1))
    allocate(q15(0:nx+1,0:ny+1,0:nz+1))
    allocate(q16(0:nx+1,0:ny+1,0:nz+1))
    allocate(q21(1:nx,1:ny,1:nz,2))
    allocate(q22(1:nx,1:ny,1:nz,2))
    allocate(q23(1:nx,1:ny,1:nz,2))
    allocate(q24(1:nx,1:ny,1:nz,2))
    allocate(q25(1:nx,1:ny,1:nz,2))
    allocate(q26(1:nx,1:ny,1:nz,2))
    allocate(dmini(nx,ny,nz))
    allocate(hxx(1:nx+1,1:nz+1))
    allocate(hyy(1:nx+1,1:nz+1))
    allocate(hzz(1:nx+1,1:nz+1))
    allocate(spa0(1:nx,1:nz))
    allocate(hxx1(1:nx+1,1:nz+1))
    allocate(hyy1(1:nx+1,1:nz+1))
    allocate(hzz1(1:nx+1,1:nz+1))
    allocate(hxx2(1:nx+1,1:nz+1))
    allocate(hyy2(1:nx+1,1:nz+1))
    allocate(hzz2(1:nx+1,1:nz+1))
    allocate(spa01(1:nx,1:nz))
    allocate(spa02(1:nx,1:nz))
    allocate(hspa(1:nx,1:nz))
    allocate(xxs(3,1:ny+1,1:nz+1,0:lb(1)+lb(2)-1)) !接收到的相邻叶排边界网格
    allocate(yys(3,1:ny+1,1:nz+1,0:lb(1)+lb(2)-1))
    allocate(zzs(3,1:ny+1,1:nz+1,0:lb(1)+lb(2)-1))
    allocate(vv(1:ny,1:nz))
    allocate(v(1:ny,1:nz,1:2,1:ny,1:nz,1:3)) !某进程网格重叠所占当前网格的面积比
    allocate(betax(ny,nz))
    allocate(betay(ny,nz))
    allocate(betaz(ny,nz))
    allocate(peb(0:ny))!出口静压
    allocate(ts1(nx,ny,nz))
    allocate(ts2(nx,ny,nz))
    allocate(ts3(nx,ny,nz))
    allocate(ts4(nx,ny,nz))
    allocate(ts5(nx,ny,nz))
    allocate(ts6(nx,ny,nz))
    allocate(s1x(nx,ny,0:nz+2)) !z方向
    allocate(s1y(nx,ny,0:nz+2))
    allocate(s1z(nx,ny,0:nz+2))
    allocate(s2x(0:nx+2,ny,nz)) !x方向
    allocate(s2y(0:nx+2,ny,nz))
    allocate(s2z(0:nx+2,ny,nz))
    allocate(s3x(nx,0:ny+2,nz)) !y方向
    allocate(s3y(nx,0:ny+2,nz))
    allocate(s3z(nx,0:ny+2,nz))
    allocate(vv0(nx,ny,nz))
    allocate(xx01(nx,ny,nz+1))
    allocate(yy01(nx,ny,nz+1))
    allocate(zz01(nx,ny,nz+1))
    allocate(xx02(nx+1,ny,nz))
    allocate(yy02(nx+1,ny,nz))
    allocate(zz02(nx+1,ny,nz))
    allocate(xx03(nx,ny+1,nz))
    allocate(yy03(nx,ny+1,nz))
    allocate(zz03(nx,ny+1,nz))
    allocate(xx0(1:nx,1:ny,1:nz))
    allocate(yy0(1:nx,1:ny,1:nz))
    allocate(zz0(1:nx,1:ny,1:nz))
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
    allocate(pvx(0:nx+1,0:ny+1,0:nz+1))
    allocate(pvy(0:nx+1,0:ny+1,0:nz+1))
    allocate(pvz(0:nx+1,0:ny+1,0:nz+1))
    allocate(vth(0:nx+1,0:ny+1,0:nz+1))
    allocate(vre(0:nx+1,0:ny+1,0:nz+1))
    allocate(p(0:nx+1,0:ny+1,0:nz+1))
    allocate(t(0:nx+1,0:ny+1,0:nz+1))
    allocate(time(nx,ny,nz))
    allocate(sri(0:nx+1,0:ny+1,0:nz+1))
    allocate(srj(0:nx+1,0:ny+1,0:nz+1))
    allocate(srk(0:nx+1,0:ny+1,0:nz+1))
    end subroutine allocation
    
subroutine ini2 !初始化流场
    use global
    implicit none
    real(8) ::vxx,vrr,vtt,vee,vth1,vre1
    integer :: num1,num2,num3
    !确定流线斜率角
    vxx=1.d0
    vtt=vxx*beta1       !展向气流角*vxx
    vrr=vxx*beta2       !周向气流角*vxx
    vee=sqrt(vxx*vxx+vrr*vrr+vtt*vtt)
    do k=1,nz
        do j=1,ny
            y1=0.25d0*(y(1,j,k)+y(1,j+1,k)+y(1,j,k+1)+y(1,j+1,k+1))
            z1=0.25d0*(z(1,j,k)+z(1,j+1,k)+z(1,j,k+1)+z(1,j+1,k+1))
            rr=sqrt(y1*y1+z1*z1)
            sir=z1/rr
            cor=y1/rr
            betax(j,k)=vxx/vee                !第一级动叶
            betay(j,k)=(vrr*cor-vtt*sir)/vee
            betaz(j,k)=(vrr*sir+vtt*cor)/vee
        end do
    end do
    rout=3.5d0*pt/ht
    dim=rout*(1.d0+0.2d0*ma(rm)*ma(rm))**(-2.5)
    q11=dim
    pp=pt*(1.d0+0.2d0*ma(rm)*ma(rm))**(-3.5)
    a=sqrt(1.4d0*pp/dim)
    en=2.5d0*pp+0.5d0*dim*ma(rm)*ma(rm)*a*a           !文件给的是绝对马赫数
    q15=en
    q16=200.*cvl0
    dimin=dim
    ppin=pp
    q12=dim*a*ma(rm)
    if(iread==0) then
        num1=mod(myid,numpp)
        if(rm==1)then
        else
            num1=num1+10
        end if
        write(id_num1,'(i3)')num1
        open(300,file='conservative_var-'//trim(adjustl(id_num1))//'myid.dat',form="unformatted")!使用定常计算中最大一个收敛背压的流场为初场
        rewind(300)
        read(300) num2      !读入继续迭代的阶段的网格层
        read(300) num3
        read(300) q11(1:nx,1:ny,1:nz)
        read(300) q12(1:nx,1:ny,1:nz)
        read(300) q13(1:nx,1:ny,1:nz)
        read(300) q14(1:nx,1:ny,1:nz)
        read(300) q15(1:nx,1:ny,1:nz)
        read(300) q16(1:nx,1:ny,1:nz)
        close(300)
        do k=1,nz   !0叶道上的q13,q14转换为当前叶道
            do j=1,ny
                do i=1,nx
                    vy=Q13(i,j,k)/Q11(i,j,k)
                    vz=Q14(i,j,k)/Q11(i,j,k)
                    y1=0.125d0*(y0(i,j,k)+y0(i+1,j,k)+y0(i,j,k+1)+y0(i+1,j,k+1)+y0(i,j+1,k)+y0(i+1,j+1,k)+y0(i,j+1,k+1)+y0(i+1,j+1,k+1))
                    z1=0.125d0*(z0(i,j,k)+z0(i+1,j,k)+z0(i,j,k+1)+z0(i+1,j,k+1)+z0(i,j+1,k)+z0(i+1,j+1,k)+z0(i,j+1,k+1)+z0(i+1,j+1,k+1))
                    rr=sqrt(y1*y1+z1*z1)
                    sir=z1/rr
                    cor=y1/rr
                    vth1=vz*cor-vy*sir    !叶高y速度  转换为0叶道上的圆柱速度
                    vre1=vz*sir+vy*cor    !周向z速度
                    
                    y1=0.125d0*(y(i,j,k)+y(i+1,j,k)+y(i,j,k+1)+y(i+1,j,k+1)+y(i,j+1,k)+y(i+1,j+1,k)+y(i,j+1,k+1)+y(i+1,j+1,k+1))!叶高y速度  转换为当前叶道上的圆柱速度
                    z1=0.125d0*(z(i,j,k)+z(i+1,j,k)+z(i,j,k+1)+z(i+1,j,k+1)+z(i,j+1,k)+z(i+1,j+1,k)+z(i,j+1,k+1)+z(i+1,j+1,k+1))
                    rr=sqrt(y1*y1+z1*z1)
                    sir=z1/rr
                    cor=y1/rr
                    vy=vre1*cor-vth1*sir
                    vz=vre1*sir+vth1*cor
                    q13(i,j,k)= q11(i,j,k)*vy
                    q14(i,j,k)= q11(i,j,k)*vz
                end do
            end do
        end do
        mm=1
        nn=2
        nitt=0
        do i=1,2
            q21(1:nx,1:ny,1:nz,i)=q11(1:nx,1:ny,1:nz)!中间量，在求导数时须要知道前两个时层值
            q22(1:nx,1:ny,1:nz,i)=q12(1:nx,1:ny,1:nz)
            q23(1:nx,1:ny,1:nz,i)=q13(1:nx,1:ny,1:nz)
            q24(1:nx,1:ny,1:nz,i)=q14(1:nx,1:ny,1:nz)
            q25(1:nx,1:ny,1:nz,i)=q15(1:nx,1:ny,1:nz)
            q26(1:nx,1:ny,1:nz,i)=q16(1:nx,1:ny,1:nz)
        end do
        if(myid==0)then
            open(4,file='convergence-inflow.dat',status="replace")
        end if
        if(rm==2 .and. myid==lb(1)*numpp+xln(0)*yln+xln(1)*yln+xln(2)-1)then
            open(7,file='convergence-outflow.dat',status="replace")
            open(8,file='pbi.dat',status="replace")
            open(11,file='eff.dat',status="replace")
        end if
        open(20,file='error-'//trim(adjustl(id_m))//'.dat',status="replace")    !如果均匀初始化，error文件重写
        if(lm==0 .and. slidm/=0 .and. yl==0)then
            open(121,file='flu'//trim(adjustl(id_rm))//'-1.dat',status="replace")
            open(122,file='flu'//trim(adjustl(id_rm))//'-2.dat',status="replace")
            open(123,file='flu'//trim(adjustl(id_rm))//'-3.dat',status="replace")
            open(124,file='flu'//trim(adjustl(id_rm))//'-4.dat',status="replace")
            open(125,file='flu'//trim(adjustl(id_rm))//'-5.dat',status="replace")
            open(126,file='flu'//trim(adjustl(id_rm))//'-6.dat',status="replace")
        end if
        if(rm==1)then!动叶，静叶前后缘每个叶道顶部都布置探针
            if(myid==3+lbm*numpp)then
                open(131,file='rl-u-'//trim(adjustl(id_m))//'.dat',status="replace")
                open(132,file='rl-p-'//trim(adjustl(id_m))//'.dat',status="replace")
            end if
            if(myid==9+lbm*numpp)then
                open(133,file='rt-u-'//trim(adjustl(id_m))//'.dat',status="replace")
                open(134,file='rt-p-'//trim(adjustl(id_m))//'.dat',status="replace")
            end if
        end if
        if(rm==2)then
            if(myid==1+lbm*numpp)then
                open(231,file='sl-u-'//trim(adjustl(id_m))//'.dat',status="replace")
                open(232,file='sl-p-'//trim(adjustl(id_m))//'.dat',status="replace")
            end if
            if(myid==8+lbm*numpp)then
                open(233,file='st-u-'//trim(adjustl(id_m))//'.dat',status="replace")
                open(234,file='st-p-'//trim(adjustl(id_m))//'.dat',status="replace")
            end if
        end if
    else if (iread==1) then
        open(myid+300,file='pause-'//trim(adjustl(id_m))//'myid.dat',form="unformatted")!unformatted使用二进制格式保存
        rewind(myid+300)
        read(myid+300) mm               !读入继续迭代的周期循环数
        read(myid+300) nn               !读入继续迭代的第n物理时步
        read(myid+300) nitt             !读入继续迭代的第nitt迭代步
        read(myid+300) q21(1:nx,1:ny,1:nz,1:2)
        read(myid+300) q22(1:nx,1:ny,1:nz,1:2)
        read(myid+300) q23(1:nx,1:ny,1:nz,1:2)
        read(myid+300) q24(1:nx,1:ny,1:nz,1:2)
        read(myid+300) q25(1:nx,1:ny,1:nz,1:2)
        read(myid+300) q26(1:nx,1:ny,1:nz,1:2)
        close(myid+300)
        q11(1:nx,1:ny,1:nz)=q21(1:nx,1:ny,1:nz,2)
        q12(1:nx,1:ny,1:nz)=q22(1:nx,1:ny,1:nz,2)
        q13(1:nx,1:ny,1:nz)=q23(1:nx,1:ny,1:nz,2)
        q14(1:nx,1:ny,1:nz)=q24(1:nx,1:ny,1:nz,2)
        q15(1:nx,1:ny,1:nz)=q25(1:nx,1:ny,1:nz,2)
        q16(1:nx,1:ny,1:nz)=q26(1:nx,1:ny,1:nz,2)
        nn=nn+1
        if(myid==0)then
            open(4,file='convergence-inflow.dat',status="old",position="append")
        end if
        if(rm==2 .and. myid==lb(1)*numpp+xln(0)*yln+xln(1)*yln+xln(2)-1)then
            open(7,file='convergence-outflow.dat',status="old",position="append")
            open(8,file='pbi.dat',status="old",position="append")
            open(11,file='eff.dat',status="old",position="append")
        end if
        open(20,file='error-'//trim(adjustl(id_m))//'.dat',status="old",position="append")
        if(lm==0 .and. slidm/=0 .and. yl==0)then
            open(121,file='flu'//trim(adjustl(id_rm))//'-1.dat',status="old",position="append")
            open(122,file='flu'//trim(adjustl(id_rm))//'-2.dat',status="old",position="append")
            open(123,file='flu'//trim(adjustl(id_rm))//'-3.dat',status="old",position="append")
            open(124,file='flu'//trim(adjustl(id_rm))//'-4.dat',status="old",position="append")
            open(125,file='flu'//trim(adjustl(id_rm))//'-5.dat',status="old",position="append")
            open(126,file='flu'//trim(adjustl(id_rm))//'-6.dat',status="old",position="append")
        end if
        if(rm==1)then!动叶，静叶前后缘每个叶道顶部都布置探针
            if(myid==3+lbm*numpp)then
                open(131,file='rl-u-'//trim(adjustl(id_m))//'.dat',status="old",position="append")
                open(132,file='rl-p-'//trim(adjustl(id_m))//'.dat',status="old",position="append")
            end if
            if(myid==9+lbm*numpp)then
                open(133,file='rt-u-'//trim(adjustl(id_m))//'.dat',status="old",position="append")
                open(134,file='rt-p-'//trim(adjustl(id_m))//'.dat',status="old",position="append")
            end if
        end if
        if(rm==2)then
            if(myid==1+lbm*numpp)then
                open(231,file='sl-u-'//trim(adjustl(id_m))//'.dat',status="old",position="append")
                open(232,file='sl-p-'//trim(adjustl(id_m))//'.dat',status="old",position="append")
            end if
            if(myid==8+lbm*numpp)then
                open(233,file='st-u-'//trim(adjustl(id_m))//'.dat',status="old",position="append")
                open(234,file='st-p-'//trim(adjustl(id_m))//'.dat',status="old",position="append")
            end if
        end if
    end if
    end subroutine ini2
    
subroutine distance   !到壁面的最短距离
    use global
    implicit none
    integer ::iil,iir,kkl,kkr,iib,iit,jjb,jjt,nx3,nx4,jtt,ny3,ny4
    real(8) ::d,dy1,dy2,dy,dz1,dz2,dz
    real(8),allocatable ::xwu(:,:),ywu(:,:),zwu(:,:),xwd(:,:),ywd(:,:),zwd(:,:),xwf(:,:),ywf(:,:),zwf(:,:),xwb(:,:),ywb(:,:),zwb(:,:) !考虑叶排的壁面坐标
    real(8),allocatable ::xd(:,:),yd(:,:),zd(:,:),xu(:,:),yu(:,:),zu(:,:),xb(:,:),yb(:,:),zb(:,:),xf(:,:),yf(:,:),zf(:,:)
    real(8),allocatable ::xx00(:,:,:),yy00(:,:,:),zz00(:,:,:)
    !*****网格体中心点坐标*****
    allocate(xx00(1:nx,1:ny,1:nz))
    allocate(yy00(1:nx,1:ny,1:nz))
    allocate(zz00(1:nx,1:ny,1:nz))
    do k=1,nz
        do j=1,ny
            do i=1,nx
                xx00(i,j,k)=0.125d0*(x(i,j,k)+x(i+1,j,k)+x(i,j,k+1)+x(i+1,j,k+1)+x(i,j+1,k)+x(i+1,j+1,k)+x(i,j+1,k+1)+x(i+1,j+1,k+1))
                yy00(i,j,k)=0.125d0*(y(i,j,k)+y(i+1,j,k)+y(i,j,k+1)+y(i+1,j,k+1)+y(i,j+1,k)+y(i+1,j+1,k)+y(i,j+1,k+1)+y(i+1,j+1,k+1))
                zz00(i,j,k)=0.125d0*(z(i,j,k)+z(i+1,j,k)+z(i,j,k+1)+z(i+1,j,k+1)+z(i,j+1,k)+z(i+1,j+1,k)+z(i,j+1,k+1)+z(i+1,j+1,k+1))
            end do
        end do
    end do
    nxyp=(it-ib+1)/xln(1)
    !*****壁面坐标*****
    !**************   y叶高方向上下壁面面中心坐标
    allocate(xwd(1:nx,1:nz))        !j=1
    allocate(ywd(1:nx,1:nz))
    allocate(zwd(1:nx,1:nz))
    allocate(xwu(1:nx,1:nz))        !j=ny
    allocate(ywu(1:nx,1:nz))
    allocate(zwu(1:nx,1:nz))
    allocate(xd(-nxyp:2*nxyp,1:nz))        !j=1
    allocate(yd(-nxyp:2*nxyp,1:nz))
    allocate(zd(-nxyp:2*nxyp,1:nz))
    allocate(xu(-nxyp:2*nxyp,1:nz))        !j=ny
    allocate(yu(-nxyp:2*nxyp,1:nz))
    allocate(zu(-nxyp:2*nxyp,1:nz))
    if(yl==0)then
        do k=1,nz
            do i=1,nx
                xwd(i,k)=0.25d0*(x(i,1,k)+x(i+1,1,k)+x(i,1,k+1)+x(i+1,1,k+1))
                ywd(i,k)=0.25d0*(y(i,1,k)+y(i+1,1,k)+y(i,1,k+1)+y(i+1,1,k+1))
                zwd(i,k)=0.25d0*(z(i,1,k)+z(i+1,1,k)+z(i,1,k+1)+z(i+1,1,k+1))
            end do
        end do
        call MPI_SENDRECV(xwd(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),760,xwu(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),860,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(ywd(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),770,ywu(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),870,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(zwd(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),780,zwu(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),880,MPI_COMM_WORLD,status,ierr)
    else if(yl==1)then
        do k=1,nz
            do i=1,nx
                xwu(i,k)=0.25d0*(x(i,ny+1,k)+x(i+1,ny+1,k)+x(i,ny+1,k+1)+x(i+1,ny+1,k+1))
                ywu(i,k)=0.25d0*(y(i,ny+1,k)+y(i+1,ny+1,k)+y(i,ny+1,k+1)+y(i+1,ny+1,k+1))
                zwu(i,k)=0.25d0*(z(i,ny+1,k)+z(i+1,ny+1,k)+z(i,ny+1,k+1)+z(i+1,ny+1,k+1))
            end do
        end do
        call MPI_SENDRECV(xwu(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),860,xwd(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),760,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(ywu(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),870,ywd(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),770,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(zwu(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),880,zwd(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),780,MPI_COMM_WORLD,status,ierr)
    end if
    if(xll>0)then
        if(xl==1)then
            j=myid-1
        else
            if(yl==0)then
                j=myid-xln(xlln-1)-1
            else
                j=myid-xln(xlln)-1
            end if
        end if
        if(rm==2 .and. xll==1)then
            nx1=24
        else
            nx1=nxyp
        end if
        call MPI_SENDRECV(xwu(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,j,70,xu(-nx1:-1,1:nz),nx1*nz,MPI_DOUBLE_PRECISION,j,60,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(ywu(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,j,71,yu(-nx1:-1,1:nz),nx1*nz,MPI_DOUBLE_PRECISION,j,61,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(zwu(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,j,72,zu(-nx1:-1,1:nz),nx1*nz,MPI_DOUBLE_PRECISION,j,62,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(xwd(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,j,73,xd(-nx1:-1,1:nz),nx1*nz,MPI_DOUBLE_PRECISION,j,63,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(ywd(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,j,74,yd(-nx1:-1,1:nz),nx1*nz,MPI_DOUBLE_PRECISION,j,64,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(zwd(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,j,75,zd(-nx1:-1,1:nz),nx1*nz,MPI_DOUBLE_PRECISION,j,65,MPI_COMM_WORLD,status,ierr)
    else
        nx1=-1
    end if
    xd(1:nx,1:nz)=xwd(1:nx,1:nz)
    yd(1:nx,1:nz)=ywd(1:nx,1:nz)
    zd(1:nx,1:nz)=zwd(1:nx,1:nz)
    xu(1:nx,1:nz)=xwu(1:nx,1:nz)
    yu(1:nx,1:nz)=ywu(1:nx,1:nz)
    zu(1:nx,1:nz)=zwu(1:nx,1:nz)
    if(xll<xln(0)+xln(1)+xln(2)-1)then
        if(xl==0)then
            if(slidm==2)then
                if(yl==0)then
                    j=myid+xln(xlln)+1
                else
                    j=myid+xln(xlln+1)+1
                end if
            else
                j=myid+1
            end if
        else
            if(yl==0)then
                j=myid+xln(xlln)+1
            else
                j=myid+xln(xlln+1)+1
            end if
        end if
        if(rm==1 .and. xll==3)then
            nx2=24
        else
            nx2=nxyp
        end if
        call MPI_SENDRECV(xwu(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,j,60,xu(nx+1:nx+nx2,1:nz),nx2*nz,MPI_DOUBLE_PRECISION,j,70,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(ywu(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,j,61,yu(nx+1:nx+nx2,1:nz),nx2*nz,MPI_DOUBLE_PRECISION,j,71,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(zwu(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,j,62,zu(nx+1:nx+nx2,1:nz),nx2*nz,MPI_DOUBLE_PRECISION,j,72,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(xwd(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,j,63,xd(nx+1:nx+nx2,1:nz),nx2*nz,MPI_DOUBLE_PRECISION,j,73,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(ywd(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,j,64,yd(nx+1:nx+nx2,1:nz),nx2*nz,MPI_DOUBLE_PRECISION,j,74,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(zwd(1:nx,1:nz),nx*nz,MPI_DOUBLE_PRECISION,j,65,zd(nx+1:nx+nx2,1:nz),nx2*nz,MPI_DOUBLE_PRECISION,j,75,MPI_COMM_WORLD,status,ierr)
    else
        nx2=0
    end if
    deallocate(xwd)
    deallocate(ywd)
    deallocate(zwd)
    deallocate(xwu)
    deallocate(ywu)
    deallocate(zwu)
    !************  z方向前后壁面面中心坐标
    if(rm==1 .and. yl==0)then
        jtt=24
    else
        jtt=ny
    end if
    allocate(xf(-nxyp:2*nxyp,-jtt:jt0))        !k=1
    allocate(yf(-nxyp:2*nxyp,-jtt:jt0))
    allocate(zf(-nxyp:2*nxyp,-jtt:jt0))
    allocate(xb(-nxyp:2*nxyp,-jtt:jt0))        !k=nz
    allocate(yb(-nxyp:2*nxyp,-jtt:jt0))
    allocate(zb(-nxyp:2*nxyp,-jtt:jt0))
    if(xlln==1)then
        allocate(xwf(1:nxyp,-jtt:jt0))        !k=1
        allocate(ywf(1:nxyp,-jtt:jt0))
        allocate(zwf(1:nxyp,-jtt:jt0))
        allocate(xwb(1:nxyp,-jtt:jt0))        !k=nz
        allocate(ywb(1:nxyp,-jtt:jt0))
        allocate(zwb(1:nxyp,-jtt:jt0))
        do j=jb,jt   !默认z方向不分块
            do i=1,nxyp
                 xwf(i,j)=0.25d0*(x(i,j,1)+x(i+1,j,1)+x(i,j+1,1)+x(i+1,j+1,1))    !前front壁面坐标1横坐标
                 ywf(i,j)=0.25d0*(y(i,j,1)+y(i+1,j,1)+y(i,j+1,1)+y(i+1,j+1,1))
                 zwf(i,j)=0.25d0*(z(i,j,1)+z(i+1,j,1)+z(i,j+1,1)+z(i+1,j+1,1))
                 xwb(i,j)=0.25d0*(x(i,j,nz+1)+x(i+1,j,nz+1)+x(i,j+1,nz+1)+x(i+1,j+1,nz+1))    !后back壁面坐标nz+1横坐标
                 ywb(i,j)=0.25d0*(y(i,j,nz+1)+y(i+1,j,nz+1)+y(i,j+1,nz+1)+y(i+1,j+1,nz+1))
                 zwb(i,j)=0.25d0*(z(i,j,nz+1)+z(i+1,j,nz+1)+z(i,j+1,nz+1)+z(i+1,j+1,nz+1))
            end do
        end do
        !**固壁区内叶片上下传递固壁值
        if(yl==0)then
            call MPI_SENDRECV(xwf(1:nxyp,jb:jt),nxyp*(jt-jb+1),MPI_DOUBLE_PRECISION,myid+xln(xlln),360,xwf(1:nxyp,jt+jb:jt+jtt),nxyp*(jtt-jb+1),MPI_DOUBLE_PRECISION,myid+xln(xlln),460,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(ywf(1:nxyp,jb:jt),nxyp*(jt-jb+1),MPI_DOUBLE_PRECISION,myid+xln(xlln),370,ywf(1:nxyp,jt+jb:jt+jtt),nxyp*(jtt-jb+1),MPI_DOUBLE_PRECISION,myid+xln(xlln),470,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(zwf(1:nxyp,jb:jt),nxyp*(jt-jb+1),MPI_DOUBLE_PRECISION,myid+xln(xlln),380,zwf(1:nxyp,jt+jb:jt+jtt),nxyp*(jtt-jb+1),MPI_DOUBLE_PRECISION,myid+xln(xlln),480,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(xwb(1:nxyp,jb:jt),nxyp*(jt-jb+1),MPI_DOUBLE_PRECISION,myid+xln(xlln),361,xwb(1:nxyp,jt+jb:jt+jtt),nxyp*(jtt-jb+1),MPI_DOUBLE_PRECISION,myid+xln(xlln),461,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(ywb(1:nxyp,jb:jt),nxyp*(jt-jb+1),MPI_DOUBLE_PRECISION,myid+xln(xlln),371,ywb(1:nxyp,jt+jb:jt+jtt),nxyp*(jtt-jb+1),MPI_DOUBLE_PRECISION,myid+xln(xlln),471,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(zwb(1:nxyp,jb:jt),nxyp*(jt-jb+1),MPI_DOUBLE_PRECISION,myid+xln(xlln),381,zwb(1:nxyp,jt+jb:jt+jtt),nxyp*(jtt-jb+1),MPI_DOUBLE_PRECISION,myid+xln(xlln),481,MPI_COMM_WORLD,status,ierr)
            ny3=jb
            ny4=jt0
        else if(yl==1)then
            call MPI_SENDRECV(xwf(1:nxyp,jb:jt),nxyp*(jt-jb+1),MPI_DOUBLE_PRECISION,myid-xln(xlln),460,xwf(1:nxyp,-jtt:-jb),nxyp*(jtt-jb+1),MPI_DOUBLE_PRECISION,myid-xln(xlln),360,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(ywf(1:nxyp,jb:jt),nxyp*(jt-jb+1),MPI_DOUBLE_PRECISION,myid-xln(xlln),470,ywf(1:nxyp,-jtt:-jb),nxyp*(jtt-jb+1),MPI_DOUBLE_PRECISION,myid-xln(xlln),370,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(zwf(1:nxyp,jb:jt),nxyp*(jt-jb+1),MPI_DOUBLE_PRECISION,myid-xln(xlln),480,zwf(1:nxyp,-jtt:-jb),nxyp*(jtt-jb+1),MPI_DOUBLE_PRECISION,myid-xln(xlln),380,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(xwb(1:nxyp,jb:jt),nxyp*(jt-jb+1),MPI_DOUBLE_PRECISION,myid-xln(xlln),461,xwb(1:nxyp,-jtt:-jb),nxyp*(jtt-jb+1),MPI_DOUBLE_PRECISION,myid-xln(xlln),361,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(ywb(1:nxyp,jb:jt),nxyp*(jt-jb+1),MPI_DOUBLE_PRECISION,myid-xln(xlln),471,ywb(1:nxyp,-jtt:-jb),nxyp*(jtt-jb+1),MPI_DOUBLE_PRECISION,myid-xln(xlln),371,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(zwb(1:nxyp,jb:jt),nxyp*(jt-jb+1),MPI_DOUBLE_PRECISION,myid-xln(xlln),481,zwb(1:nxyp,-jtt:-jb),nxyp*(jtt-jb+1),MPI_DOUBLE_PRECISION,myid-xln(xlln),381,MPI_COMM_WORLD,status,ierr)
            ny3=-jtt
            ny4=jt
        end if
        xb(1:nxyp,ny3:ny4)=xwb(1:nxyp,ny3:ny4)
        yb(1:nxyp,ny3:ny4)=ywb(1:nxyp,ny3:ny4)
        zb(1:nxyp,ny3:ny4)=zwb(1:nxyp,ny3:ny4)
        xf(1:nxyp,ny3:ny4)=xwf(1:nxyp,ny3:ny4)
        yf(1:nxyp,ny3:ny4)=ywf(1:nxyp,ny3:ny4)
        zf(1:nxyp,ny3:ny4)=zwf(1:nxyp,ny3:ny4)
        !**向左右流场和固壁区传递固壁值
        if(xl==0)then
            if(yl==0)then
                j=0
            else
                j=1
            end if
            do i=1,xln(0)
                call MPI_SEND(ny3,1,MPI_INTEGER,myid-xln(j)-i,210,MPI_COMM_WORLD,ierr)
                call MPI_SEND(ny4,1,MPI_INTEGER,myid-xln(j)-i,220,MPI_COMM_WORLD,ierr)
                call MPI_SEND(xwb(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-xln(j)-i,230,MPI_COMM_WORLD,ierr)
                call MPI_SEND(ywb(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-xln(j)-i,240,MPI_COMM_WORLD,ierr)
                call MPI_SEND(zwb(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-xln(j)-i,250,MPI_COMM_WORLD,ierr)
                call MPI_SEND(xwf(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-xln(j)-i,260,MPI_COMM_WORLD,ierr)
                call MPI_SEND(ywf(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-xln(j)-i,270,MPI_COMM_WORLD,ierr)
                call MPI_SEND(zwf(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-xln(j)-i,280,MPI_COMM_WORLD,ierr)
            end do
            call MPI_SENDRECV(xwb(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+1,231,xb(nxyp+1:2*nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+1,241,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(ywb(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+1,232,yb(nxyp+1:2*nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+1,242,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(zwb(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+1,233,zb(nxyp+1:2*nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+1,243,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(xwf(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+1,234,xf(nxyp+1:2*nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+1,244,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(ywf(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+1,235,yf(nxyp+1:2*nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+1,245,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(zwf(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+1,236,zf(nxyp+1:2*nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+1,246,MPI_COMM_WORLD,status,ierr)
            nx3=1
            nx4=2*nxyp
        end if
        if(xl==xln(1)-1)then
            if(yl==0)then
                j=1
            else
                j=2
            end if
            do i=1,xln(2)
                call MPI_SEND(ny3,1,MPI_INTEGER,myid+xln(j)+i,210,MPI_COMM_WORLD,ierr)
                call MPI_SEND(ny4,1,MPI_INTEGER,myid+xln(j)+i,220,MPI_COMM_WORLD,ierr)
                call MPI_SEND(xwb(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+xln(j)+i,230,MPI_COMM_WORLD,ierr)
                call MPI_SEND(ywb(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+xln(j)+i,240,MPI_COMM_WORLD,ierr)
                call MPI_SEND(zwb(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+xln(j)+i,250,MPI_COMM_WORLD,ierr)
                call MPI_SEND(xwf(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+xln(j)+i,260,MPI_COMM_WORLD,ierr)
                call MPI_SEND(ywf(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+xln(j)+i,270,MPI_COMM_WORLD,ierr)
                call MPI_SEND(zwf(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid+xln(j)+i,280,MPI_COMM_WORLD,ierr)
            end do
            call MPI_SENDRECV(xwb(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-1,241,xb(-nxyp:-1,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-1,231,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(ywb(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-1,242,yb(-nxyp:-1,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-1,232,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(zwb(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-1,243,zb(-nxyp:-1,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-1,233,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(xwf(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-1,244,xf(-nxyp:-1,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-1,234,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(ywf(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-1,245,yf(-nxyp:-1,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-1,235,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(zwf(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-1,246,zf(-nxyp:-1,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,myid-1,236,MPI_COMM_WORLD,status,ierr)
            nx3=-nxyp
            nx4=nxyp
        end if
        deallocate(xwf)
        deallocate(ywf)
        deallocate(zwf)
        deallocate(xwb)
        deallocate(ywb)
        deallocate(zwb)
    else
        if(xlln==0)then
            if(yl==0)then
                j=xln(0)*yln+lbm*numpp
            else if(yl==1)then
                j=xln(0)*yln+xln(1)+lbm*numpp
            end if
        else if(xlln==2)then
            if(yl==0)then
                j=xln(0)*yln+xln(1)-1+lbm*numpp
            else if(yl==1)then
                j=xln(0)*yln+xln(1)*yln-1+lbm*numpp
            end if
        end if
        call MPI_RECV(ny3,1,MPI_INTEGER,j,210,MPI_COMM_WORLD,status,ierr)
        call MPI_RECV(ny4,1,MPI_INTEGER,j,220,MPI_COMM_WORLD,status,ierr)
        call MPI_RECV(xb(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,j,230,MPI_COMM_WORLD,status,ierr)
        call MPI_RECV(yb(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,j,240,MPI_COMM_WORLD,status,ierr)
        call MPI_RECV(zb(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,j,250,MPI_COMM_WORLD,status,ierr)
        call MPI_RECV(xf(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,j,260,MPI_COMM_WORLD,status,ierr)
        call MPI_RECV(yf(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,j,270,MPI_COMM_WORLD,status,ierr)
        call MPI_RECV(zf(1:nxyp,ny3:ny4),nxyp*(ny4-ny3+1),MPI_DOUBLE_PRECISION,j,280,MPI_COMM_WORLD,status,ierr)
    end if
    !*****求到壁面最短距离*****
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    dy=1.d20
                    !求距离上下固壁距离,先确定某空间点下对应的固壁点坐标范围
                    iil=i-10
                    iil=max(iil,-nx1)
                    iir=iil+20
                    iir=min(iir,nx+nx2)

                    kkl=k-10
                    kkl=max(kkl,1)
                    kkr=kkl+20
                    kkr=min(kkr,nz)
                    
                    do kk=kkl,kkr
                        do ii=iil,iir
                            if(ii==0)then
                            else
                                dy1=(xx00(i,j,k)-xd(ii,kk))**2+(yy00(i,j,k)-yd(ii,kk))**2+(zz00(i,j,k)-zd(ii,kk))**2
                                dy2=(xx00(i,j,k)-xu(ii,kk))**2+(yy00(i,j,k)-yu(ii,kk))**2+(zz00(i,j,k)-zu(ii,kk))**2
                                d=min(dy1,dy2)
                                if(d<dy)then
                                    dy=d
                                end if
                            end if
                        end do
                    end do
                    !求距离z方向的前后固壁距离,同上方法，只是要讨论下
                    dz=1.d20
                    if(xlln==1)then
                        iib=i-10
                        iib=max(iib,nx3)
                        iit=iib+20
                        iit=min(iit,nx4)
                    else if(xlln==0)then
                        iib=1
                        iit=iib+20
                    else if(xlln==2)then
                        iib=nxyp-20
                        iit=nxyp
                    end if
                    if(j<jb) then
                        jjb=jb
                        jjt=jb+20
                    else if (j>=jb .and. j<=jt) then
                        jjb=j-10
                        jjb=max(jjb,ny3)
                        jjt=jjb+20
                        jjt=min(jjt,ny4)
                    else
                        jjb=jt-20
                        jjt=jt
                    end if
                    
                    do jj=jjb,jjt
                        do ii=iib,iit
                            if(ii==0 .or. jj==0)then
                            else
                                dz1=(xx00(i,j,k)-xf(ii,jj))**2+(yy00(i,j,k)-yf(ii,jj))**2+(zz00(i,j,k)-zf(ii,jj))**2
                                dz2=(xx00(i,j,k)-xb(ii,jj))**2+(yy00(i,j,k)-yb(ii,jj))**2+(zz00(i,j,k)-zb(ii,jj))**2
                                d=min(dz1,dz2)
                                if (d<dz) then
                                    dz=d
                                end if
                            end if
                        end do
                    end do
                    dmini(i,j,k)=sqrt(min(dy,dz)) !dy,dz求较小值
                end do
            end do
        end do
        deallocate(xx00)
        deallocate(yy00)
        deallocate(zz00)
        deallocate(xd)
        deallocate(yd)
        deallocate(zd)
        deallocate(xu)
        deallocate(yu)
        deallocate(zu)
        deallocate(xf)
        deallocate(yf)
        deallocate(zf)
        deallocate(xb)
        deallocate(yb)
        deallocate(zb)
    end subroutine distance
    
subroutine geo
    use global
    implicit none
    real(8) :: temp,x24,y24,z24,x31,y31,z31
    
    do k=1,nz+1
        do j=1,ny
            do i=1,nx
               x24=x(i,j+1,k)-x(i+1,j,k)
               y24=y(i,j+1,k)-y(i+1,j,k)
               z24=z(i,j+1,k)-z(i+1,j,k)

               x31=x(i+1,j+1,k)-x(i,j,k)
               y31=y(i+1,j+1,k)-y(i,j,k)
               z31=z(i+1,j+1,k)-z(i,j,k)

               s1x(i,j,k)=0.5d0*(y24*z31-z24*y31)   !前后面z方向1
               s1y(i,j,k)=0.5d0*(z24*x31-x24*z31)
               s1z(i,j,k)=0.5d0*(x24*y31-y24*x31)
               
               xx01(i,j,k)=0.25d0*(x(i,j,k)+x(i,j+1,k)+x(i+1,j,k)+x(i+1,j+1,k))
               yy01(i,j,k)=0.25d0*(y(i,j,k)+y(i,j+1,k)+y(i+1,j,k)+y(i+1,j+1,k))
               zz01(i,j,k)=0.25d0*(z(i,j,k)+z(i,j+1,k)+z(i+1,j,k)+z(i+1,j+1,k))
            end do
        end do
    end do
    do k=1,nz
        do j=1,ny
            do i=1,nx+1
               x24=x(i,j+1,k+1)-x(i,j,k)
               y24=y(i,j+1,k+1)-y(i,j,k)
               z24=z(i,j+1,k+1)-z(i,j,k)

               x31=x(i,j+1,k)-x(i,j,k+1)
               y31=y(i,j+1,k)-y(i,j,k+1)
               z31=z(i,j+1,k)-z(i,j,k+1)

               s2x(i,j,k)=0.5d0*(y24*z31-z24*y31)  !左右面x方向2
               s2y(i,j,k)=0.5d0*(z24*x31-x24*z31)
               s2z(i,j,k)=0.5d0*(x24*y31-y24*x31)
               
               xx02(i,j,k)=0.25d0*(x(i,j,k)+x(i,j+1,k)+x(i,j,k+1)+x(i,j+1,k+1))
               yy02(i,j,k)=0.25d0*(y(i,j,k)+y(i,j+1,k)+y(i,j,k+1)+y(i,j+1,k+1))
               zz02(i,j,k)=0.25d0*(z(i,j,k)+z(i,j+1,k)+z(i,j,k+1)+z(i,j+1,k+1))
            end do
        end do
    end do
    do k=1,nz
        do j=1,ny+1
            do i=1,nx
               x24=x(i+1,j,k+1)-x(i,j,k)
               y24=y(i+1,j,k+1)-y(i,j,k)
               z24=z(i+1,j,k+1)-z(i,j,k)

               x31=x(i,j,k+1)-x(i+1,j,k)
               y31=y(i,j,k+1)-y(i+1,j,k)
               z31=z(i,j,k+1)-z(i+1,j,k)

               s3x(i,j,k)=0.5d0*(y24*z31-z24*y31)  !上下面y方向3
               s3y(i,j,k)=0.5d0*(z24*x31-x24*z31)
               s3z(i,j,k)=0.5d0*(x24*y31-y24*x31)
               
               xx03(i,j,k)=0.25d0*(x(i,j,k)+x(i+1,j,k)+x(i,j,k+1)+x(i+1,j,k+1))
               yy03(i,j,k)=0.25d0*(y(i,j,k)+y(i+1,j,k)+y(i,j,k+1)+y(i+1,j,k+1))
               zz03(i,j,k)=0.25d0*(z(i,j,k)+z(i+1,j,k)+z(i,j,k+1)+z(i+1,j,k+1))
            end do
        end do
    end do
    do k=1,nz
        do j=1,ny
            do i=1,nx
               x24=x(i+1,j+1,k+1)-x(i,j,k)
               y24=y(i+1,j+1,k+1)-y(i,j,k)
               z24=z(i+1,j+1,k+1)-z(i,j,k)

               vv0(i,j,k)=-(x24*(s1x(i,j,k)+s2x(i,j,k)+s3x(i,j,k))  &
                           +y24*(s1y(i,j,k)+s2y(i,j,k)+s3y(i,j,k))  &
                           +z24*(s1z(i,j,k)+s2z(i,j,k)+s3z(i,j,k)))/3.d0
               if(vv0(i,j,k)<=0.d0) then
                    write(*,*) 'vv0<=0',i,j,k
               end if
               xx0(i,j,k)=0.125d0*(x(i,j,k)+x(i+1,j,k)+x(i,j,k+1)+x(i+1,j,k+1)+x(i,j+1,k)+x(i+1,j+1,k)+x(i,j+1,k+1)+x(i+1,j+1,k+1))
               yy0(i,j,k)=0.125d0*(y(i,j,k)+y(i+1,j,k)+y(i,j,k+1)+y(i+1,j,k+1)+y(i,j+1,k)+y(i+1,j+1,k)+y(i,j+1,k+1)+y(i+1,j+1,k+1))
               zz0(i,j,k)=0.125d0*(z(i,j,k)+z(i+1,j,k)+z(i,j,k+1)+z(i+1,j,k+1)+z(i,j+1,k)+z(i+1,j+1,k)+z(i,j+1,k+1)+z(i+1,j+1,k+1))
            end do
        end do
    end do
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
    !虚拟网格顶点坐标值,及交界面第一层网格体积
    if(slidm==1)then
        x(nx+2,1:ny+1,1:nz+1)=2.*x(nx+1,1:ny+1,1:nz+1)-x(nx,1:ny+1,1:nz+1) !网格顶点坐标
        y(nx+2,1:ny+1,1:nz+1)=2.*y(nx+1,1:ny+1,1:nz+1)-y(nx,1:ny+1,1:nz+1)
        z(nx+2,1:ny+1,1:nz+1)=2.*z(nx+1,1:ny+1,1:nz+1)-z(nx,1:ny+1,1:nz+1)
        vv(1:ny,1:nz)=vv0(nx,1:ny,1:nz)
    else if(slidm==2)then
        x(0,1:ny+1,1:nz+1)=2.*x(1,1:ny+1,1:nz+1)-x(2,1:ny+1,1:nz+1)
        y(0,1:ny+1,1:nz+1)=2.*y(1,1:ny+1,1:nz+1)-y(2,1:ny+1,1:nz+1)
        z(0,1:ny+1,1:nz+1)=2.*z(1,1:ny+1,1:nz+1)-z(2,1:ny+1,1:nz+1)
        vv(1:ny,1:nz)=vv0(1,1:ny,1:nz)
    end if
    end subroutine geo
    
subroutine spaa  !由于y方向分区，为后处理输出打基础,求出当前进程下叶高百分比对应的实际距离值
    use global
    implicit none
    real(8),allocatable ::xff(:,:,:),yff(:,:,:),zff(:,:,:),sss(:),hx(:),hy(:),hz(:),hr(:)
    real(8) :: t1,t2,temp
    
    ny0=ny*yln
    allocate(xff(1:nx+1,1:ny0+1,1:nz+1))!分配得yln=1下的网格
    allocate(yff(1:nx+1,1:ny0+1,1:nz+1))
    allocate(zff(1:nx+1,1:ny0+1,1:nz+1))
    if(rm==1)then
        nx1=nxyp!当前进程前的每进程网格数
        xff(1:nx+1,1:ny0+1,1:nz+1)=xf0(xll*nx1+1:xll*nx1+nx+1,1:ny0+1,zl*nz+1:(zl+1)*nz+1)
        yff(1:nx+1,1:ny0+1,1:nz+1)=yf0(xll*nx1+1:xll*nx1+nx+1,1:ny0+1,zl*nz+1:(zl+1)*nz+1)
        zff(1:nx+1,1:ny0+1,1:nz+1)=zf0(xll*nx1+1:xll*nx1+nx+1,1:ny0+1,zl*nz+1:(zl+1)*nz+1)
    else
        nx1=24!当前进程前的每进程网格数
        xff(1:nx+1,1:ny0+1,1:nz+1)=xf0((xll-1)*nx+nx1+1:xll*nx+nx1+1,1:ny0+1,zl*nz+1:(zl+1)*nz+1)
        yff(1:nx+1,1:ny0+1,1:nz+1)=yf0((xll-1)*nx+nx1+1:xll*nx+nx1+1,1:ny0+1,zl*nz+1:(zl+1)*nz+1)
        zff(1:nx+1,1:ny0+1,1:nz+1)=zf0((xll-1)*nx+nx1+1:xll*nx+nx1+1,1:ny0+1,zl*nz+1:(zl+1)*nz+1)
    end if
    temp=2.d0/dble(lb(rm))*pi*lm   !旋转得整周
    cor=cos(temp)                                       !0             ^
    sir=sin(temp)                                       !1             ^        叶片往上旋转
    do k=1,nz+1                                         !nz            ^
        do j=1,ny0+1                                      !nz+1          1
           do i=1,nx+1
               t1=yff(i,j,k)
               t2=zff(i,j,k)
               yff(i,j,k)=t1*cor+t2*sir !按叶片旋转方向旋转,若相反则是为了求0、nz+1坐标，1-nz的顺序与叶片旋转方向相反
               zff(i,j,k)=-t1*sir+t2*cor
           end do
       end do
    end do
    allocate(sss(1:ny0+1))
    allocate(hx(1:ny0+1))
    allocate(hy(1:ny0+1))
    allocate(hz(1:ny0+1))
    allocate(hr(1:ny0+1))
    if(spa1<50.)then!只是个粗略判断标准
        l=0
    else
        l=1
    end if
    if(spa2<50.)then!只是个粗略判断标准
        kk=0
    else
        kk=1
    end if
    do k=1,nz+1
        do i=1,nx+1
            sss(1)=0.d0
            do j=1,ny0+1
                hx(j)=xff(i,j,k)
                hy(j)=yff(i,j,k)
                hz(j)=zff(i,j,k)
                hr(j)=sqrt(hy(j)*hy(j)+hz(j)*hz(j))
                if(j>1)then
                    sss(j)=sss(j-1)+sqrt((hx(j)-hx(j-1))**2+(hr(j)-hr(j-1))**2) !以上为计算坐标各网格中心坐标
                end if
            end do
            if(l==0 .and. yl==l)then
                t1=sss(ny0+1)*spa1/100.
                call wl(sss(1:ny0+1),hx(1:ny0+1),t1,hxx1(i,k),ny0+1,1)
                call wl(sss(1:ny0+1),hy(1:ny0+1),t1,hyy1(i,k),ny0+1,1)
                call wl(sss(1:ny0+1),hz(1:ny0+1),t1,hzz1(i,k),ny0+1,1)
            else if (l==1 .and. yl==l)then
                t1=sss(ny0+1)*spa1/100.
                call wl(sss(1:ny0+1),hx(1:ny0+1),t1,hxx1(i,k),ny0+1,1)
                call wl(sss(1:ny0+1),hy(1:ny0+1),t1,hyy1(i,k),ny0+1,1)
                call wl(sss(1:ny0+1),hz(1:ny0+1),t1,hzz1(i,k),ny0+1,1)
            end if
            if(kk==0 .and. yl==kk)then
                t2=sss(ny0+1)*spa2/100.
                call wl(sss(1:ny0+1),hx(1:ny0+1),t2,hxx2(i,k),ny0+1,1)
                call wl(sss(1:ny0+1),hy(1:ny0+1),t2,hyy2(i,k),ny0+1,1)
                call wl(sss(1:ny0+1),hz(1:ny0+1),t2,hzz2(i,k),ny0+1,1)
            else if(kk==1 .and. yl==kk)then
                t2=sss(ny0+1)*spa2/100.
                call wl(sss(1:ny0+1),hx(1:ny0+1),t2,hxx2(i,k),ny0+1,1)
                call wl(sss(1:ny0+1),hy(1:ny0+1),t2,hyy2(i,k),ny0+1,1)
                call wl(sss(1:ny0+1),hz(1:ny0+1),t2,hzz2(i,k),ny0+1,1)
            end if
        end do
    end do
    do k=1,nz
        do i=1,nx
            sss(1)=0.d0
            do j=1,ny0
                hx(j)=0.125d0*(xff(i,j,k)+xff(i+1,j,k)+xff(i,j,k+1)+xff(i+1,j,k+1)+xff(i,j+1,k)+xff(i+1,j+1,k)+xff(i,j+1,k+1)+xff(i+1,j+1,k+1))
                y1=0.125d0*(yff(i,j,k)+yff(i+1,j,k)+yff(i,j,k+1)+yff(i+1,j,k+1)+yff(i,j+1,k)+yff(i+1,j+1,k)+yff(i,j+1,k+1)+yff(i+1,j+1,k+1))
                z1=0.125d0*(zff(i,j,k)+zff(i+1,j,k)+zff(i,j,k+1)+zff(i+1,j,k+1)+zff(i,j+1,k)+zff(i+1,j+1,k)+zff(i,j+1,k+1)+zff(i+1,j+1,k+1))
                hr(j)=sqrt(y1*y1+z1*z1)
                if(j>1)then
                    sss(j)=sss(j-1)+sqrt((hx(j)-hx(j-1))**2+(hr(j)-hr(j-1))**2) !以上为计算坐标各网格中心坐标
                end if
            end do
            hspa(i,k)=sss(ny+1)
            if(l==0 .and. yl==l)then
                spa01(i,k)=sss(ny0)*spa1/100.
            else if (l==1 .and. yl==l)then
                spa01(i,k)=sss(ny0)*spa1/100.
            end if
            if(kk==0 .and. yl==kk)then
                spa02(i,k)=sss(ny0)*spa2/100.
            else if (kk==1 .and. yl==kk)then
                spa02(i,k)=sss(ny0)*spa2/100.
            end if
        end do
    end do
    !为求物理量沿叶高分布打基础
    if(rm==2 .and. xll==4 .and. lm==0)then
        allocate(temp0(1:ny0))
        i=nx
        k=1
        do j=1,ny0
            y1=0.125d0*(yff(i,j,k)+yff(i+1,j,k)+yff(i,j,k+1)+yff(i+1,j,k+1)+yff(i,j+1,k)+yff(i+1,j+1,k)+yff(i,j+1,k+1)+yff(i+1,j+1,k+1))
            z1=0.125d0*(zff(i,j,k)+zff(i+1,j,k)+zff(i,j,k+1)+zff(i+1,j,k+1)+zff(i,j+1,k)+zff(i+1,j+1,k)+zff(i,j+1,k+1)+zff(i+1,j+1,k+1))
            hr(j)=sqrt(y1*y1+z1*z1)
        end do
        do j=1,ny0
            temp0(j)=(hr(j)-hr(1))/(hr(ny0)-hr(1))*100.d0
        end do
    end if
    deallocate(xf0)
    deallocate(yf0)
    deallocate(zf0)
    deallocate(xff)
    deallocate(yff)
    deallocate(zff)
    deallocate(sss)
    deallocate(hx)
    deallocate(hy)
    deallocate(hz)
    deallocate(hr)
    end subroutine spaa
    
    subroutine slid  !各进程将其交接面真实顶点（2层）经坐标转换后，传递给相邻叶排，与相邻叶排的虚拟网格裁剪
    use global
    implicit none
    real(8) ::temp
    real(8),allocatable :: ys(:,:,:),zs(:,:,:)
  
    temp=delt*dble(n-1+(m-1)*nt)*rpm(1)!当前外、中循环下旋转坐标转过的角度。rpm为负值，故后面动叶要取负号
    if(rm==1)then
        i=nx-1
        temp=-temp  !上游动叶的网格顶点映射到静叶上，所以要相对坐标->绝对坐标a,按叶片旋转方向转
    else if(rm==2)then
        i=1
        temp=temp   !下游静叶的网格顶点映射到动叶上，所以要绝对坐标->相绝对坐标r,按叶片旋转方向反转
    end if
    allocate(ys(1:nx+1,1:ny+1,1:nz+1))
    allocate(zs(1:nx+1,1:ny+1,1:nz+1))
    cor=cos(temp)                           !n从2开始，故刚开始计算时已经移动了0.05个通道
    sir=sin(temp)
    do l=i,i+2
        do k=1,nz+1
            do j=1,ny+1
                ys(l,j,k)=y(l,j,k)*cor+z(l,j,k)*sir
                zs(l,j,k)=-y(l,j,k)*sir+z(l,j,k)*cor
            end do
        end do
    end do
    xxs=0.
    yys=0.
    zzs=0.
    if(rm==1)then
        do j=yl+lb(rm)*numpp,yl+lb(rm)*numpp+(lb(rm+1)-1)*numpp,numpp
            k=j/numpp
            call MPI_SENDRECV(x(i:i+2,1:ny+1,1:nz+1),3*(ny+1)*(nz+1),MPI_DOUBLE_PRECISION,j,5,xxs(1:3,1:ny+1,1:nz+1,k),3*(ny+1)*(nz+1),MPI_DOUBLE_PRECISION,j,15,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(ys(i:i+2,1:ny+1,1:nz+1),3*(ny+1)*(nz+1),MPI_DOUBLE_PRECISION,j,6,yys(1:3,1:ny+1,1:nz+1,k),3*(ny+1)*(nz+1),MPI_DOUBLE_PRECISION,j,16,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(zs(i:i+2,1:ny+1,1:nz+1),3*(ny+1)*(nz+1),MPI_DOUBLE_PRECISION,j,7,zzs(1:3,1:ny+1,1:nz+1,k),3*(ny+1)*(nz+1),MPI_DOUBLE_PRECISION,j,17,MPI_COMM_WORLD,status,ierr)
        end do
    else if(rm==2)then
        do j=yl+2*xln(1)*yln,yl+2*xln(1)*yln+(lb(rm-1)-1)*numpp,numpp
            k=j/numpp
            call MPI_SENDRECV(x(i:i+2,1:ny+1,1:nz+1),3*(ny+1)*(nz+1),MPI_DOUBLE_PRECISION,j,15,xxs(1:3,1:ny+1,1:nz+1,k),3*(ny+1)*(nz+1),MPI_DOUBLE_PRECISION,j,5,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(ys(i:i+2,1:ny+1,1:nz+1),3*(ny+1)*(nz+1),MPI_DOUBLE_PRECISION,j,16,yys(1:3,1:ny+1,1:nz+1,k),3*(ny+1)*(nz+1),MPI_DOUBLE_PRECISION,j,6,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(zs(i:i+2,1:ny+1,1:nz+1),3*(ny+1)*(nz+1),MPI_DOUBLE_PRECISION,j,17,zzs(1:3,1:ny+1,1:nz+1,k),3*(ny+1)*(nz+1),MPI_DOUBLE_PRECISION,j,7,MPI_COMM_WORLD,status,ierr)
        end do
    end if
    deallocate(ys)
    deallocate(zs)
    end subroutine slid
    
    subroutine overlapzone  !当前进程虚拟网格被相邻叶排真实网格（2层）裁剪，找到与其重叠的是相邻叶排的哪几个进程,并找到重叠的是哪个k部分
    use global
    implicit none
    
    if(rm==1)then
        i=nx+1
        myidl=yl+lb(rm)*numpp
        myidr=yl+lb(rm)*numpp+(lb(rm+1)-1)*numpp
    else if(rm==2)then
        i=0
        myidl=yl+2*xln(1)*yln
        myidr=yl+2*xln(1)*yln+(lb(rm-1)-1)*numpp
    end if
    ssum=0
    cli=0
    do l=myidl,myidr,numpp
        k=l/numpp
        j=0
        call clip(xxs(1,1,1,k),xxs(1,ny+1,1,k),xxs(1,ny+1,nz+1,k),xxs(1,1,nz+1,k),xxs(3,1,1,k),xxs(3,ny+1,1,k),xxs(3,ny+1,nz+1,k),xxs(3,1,nz+1,k),&
            yys(1,1,1,k),yys(1,ny+1,1,k),yys(1,ny+1,nz+1,k),yys(1,1,nz+1,k),yys(3,1,1,k),yys(3,ny+1,1,k),yys(3,ny+1,nz+1,k),yys(3,1,nz+1,k),&
            zzs(1,1,1,k),zzs(1,ny+1,1,k),zzs(1,ny+1,nz+1,k),zzs(1,1,nz+1,k),zzs(3,1,1,k),zzs(3,ny+1,1,k),zzs(3,ny+1,nz+1,k),zzs(3,1,nz+1,k),&
            x(i,1,1),x(i,ny+1,1),x(i,ny+1,nz+1),x(i,1,nz+1),x(i+1,1,1),x(i+1,ny+1,1),x(i+1,ny+1,nz+1),x(i+1,1,nz+1),&
            y(i,1,1),y(i,ny+1,1),y(i,ny+1,nz+1),y(i,1,nz+1),y(i+1,1,1),y(i+1,ny+1,1),y(i+1,ny+1,nz+1),y(i+1,1,nz+1),&
            z(i,1,1),z(i,ny+1,1),z(i,ny+1,nz+1),z(i,1,nz+1),z(i+1,1,1),z(i+1,ny+1,1),z(i+1,ny+1,nz+1),z(i+1,1,nz+1),j)!判断两进程是否重叠，输出j
        if(j==1)then
            ssum=ssum+1 !重叠个数
            cli(ssum)=l !重叠的进程号，从小到大排列
        end if
        if(ssum/=0 .and. cli(1)/=myidl .and. j==0)then!因为Myidl与myidr会连在一起，故若myidl重叠，后面不重叠，还是不能退出，还要判断myidr。
            exit                                     !但若Myidl不重叠，后面的第2，或更多进程才重叠，则再不重叠就不用再判断了
        end if
    end do
    !cli按进程号由小到大输出，一般与小进程号（myidl）的重叠部分对应<当前进程>的大k,中进程号（myidm）对应中k,大进程号(myidr)对应小k。除了myidl和myidr外。
    if(rm==1)then  !动叶占有2~3个静叶
        if(ssum==2)then
            myidl=cli(1)
            myidr=cli(2)
            if(cli(1)==yl+lb(rm)*numpp .and. cli(2)==yl+lb(rm)*numpp+(lb(rm+1)-1)*numpp)then
                myidl=cli(2)
                myidr=cli(1)
            end if
        else if(ssum==3)then
            myidl=cli(1)
            myidm=cli(2)
            myidr=cli(3)
            if(cli(1)==yl+lb(rm)*numpp .and. cli(2)==yl+(lb(rm)+1)*numpp .and. cli(3)==yl+lb(rm)*numpp+(lb(rm+1)-1)*numpp)then
                myidl=cli(3)
                myidm=cli(1)
                myidr=cli(2)
            else if(cli(1)==yl+lb(rm)*numpp .and. cli(2)==yl+lb(rm)*numpp+(lb(rm+1)-2)*numpp .and. cli(3)==yl+lb(rm)*numpp+(lb(rm+1)-1)*numpp)then
                myidl=cli(2)
                myidm=cli(3)
                myidr=cli(1)
            end if
        end if
    else if(rm==2)then  !静叶占有1~2个静叶
        if(ssum==1)then
            myidm=cli(1)
        else if(ssum==2)then
            myidl=cli(1)
            myidr=cli(2)
            if(cli(1)==yl+2*xln(1)*yln .and. cli(2)==yl+2*xln(1)*yln+(lb(rm-1)-1)*numpp)then
                myidl=cli(2)
                myidr=cli(1)
            end if
        end if
    end if
    end subroutine overlapzone
    
   subroutine overlapgrid  !当前进程每个虚拟网格被相邻叶排进程每个真实网格拆分的面积比
    use global
    implicit none  !进程在j边缘加密,若相邻进程此位置不加密,则出现圆弧对直线的问题,故不能只对同一j,扩大范围到j+-1
    real(8) :: overv,overvv,av
    integer :: num,num1,ll,lll,iii,jjj,sum,xk,yk,zk,xj,yj
    integer :: xk1(2),kkk(2)
    
    v=0.
    do l=1,ssum
        num=cli(l)
        num1=num/numpp
        if(num==myidl)then
            xk=nz
            yk=1
            zk=-1
        else if(num==myidr .or. num==myidm)then
            xk=1
            yk=nz
            zk=1
        end if
        do j=1,ny
            if(j==1)then
                xj=1
                yj=2
            else if(j==ny)then
                xj=ny-1
                yj=ny
            else
                xj=j-1
                yj=j+1
            end if
            if(num==myidl)then
                xk1(1)=nz
                xk1(2)=nz
            else if(num==myidr .or. num==myidm)then
                xk1(1)=1
                xk1(2)=1
            end if
            jj=j
            iii=0
            do k=xk,yk,zk
                kkk=0
                do ii=1,2
                    lll=0
                    do kk=xk1(ii),yk,zk
                        overv=0.
                        call clipv(xxs(ii,jj,kk,num1),xxs(ii,jj+1,kk,num1),xxs(ii,jj+1,kk+1,num1),xxs(ii,jj,kk+1,num1),xxs(ii+1,jj,kk,num1),xxs(ii+1,jj+1,kk,num1),xxs(ii+1,jj+1,kk+1,num1),xxs(ii+1,jj,kk+1,num1),&
                            yys(ii,jj,kk,num1),yys(ii,jj+1,kk,num1),yys(ii,jj+1,kk+1,num1),yys(ii,jj,kk+1,num1),yys(ii+1,jj,kk,num1),yys(ii+1,jj+1,kk,num1),yys(ii+1,jj+1,kk+1,num1),yys(ii+1,jj,kk+1,num1),&
                            zzs(ii,jj,kk,num1),zzs(ii,jj+1,kk,num1),zzs(ii,jj+1,kk+1,num1),zzs(ii,jj,kk+1,num1),zzs(ii+1,jj,kk,num1),zzs(ii+1,jj+1,kk,num1),zzs(ii+1,jj+1,kk+1,num1),zzs(ii+1,jj,kk+1,num1),&
                            x(i,j,k),x(i,j+1,k),x(i,j+1,k+1),x(i,j,k+1),x(i+1,j,k),x(i+1,j+1,k),x(i+1,j+1,k+1),x(i+1,j,k+1),&
                            y(i,j,k),y(i,j+1,k),y(i,j+1,k+1),y(i,j,k+1),y(i+1,j,k),y(i+1,j+1,k),y(i+1,j+1,k+1),y(i+1,j,k+1),&
                            z(i,j,k),z(i,j+1,k),z(i,j+1,k+1),z(i,j,k+1),z(i+1,j,k),z(i+1,j+1,k),z(i+1,j+1,k+1),z(i+1,j,k+1),overv)
                        if(overv/=0.)then
                            iii=1   !当前进程有与相邻叶排进程重叠过
                            lll=lll+1
                            if(lll==1)then
                                xk1(ii)=kk!取第1次重叠时的kk值，决定do ii循环里下一次kk的起始值，当前do kk循环不影响，决定了相邻叶排进程被裁剪循环的起始点位置。
                            end if
                            v(j,k,ii,jj,kk,l)=overv/vv(j,k)
                            if(xj/=j)then
                                jjj=xj
                                overvv=0.
                                call clipv(xxs(ii,jjj,kk,num1),xxs(ii,jjj+1,kk,num1),xxs(ii,jjj+1,kk+1,num1),xxs(ii,jjj,kk+1,num1),xxs(ii+1,jjj,kk,num1),xxs(ii+1,jjj+1,kk,num1),xxs(ii+1,jjj+1,kk+1,num1),xxs(ii+1,jjj,kk+1,num1),&
                                    yys(ii,jjj,kk,num1),yys(ii,jjj+1,kk,num1),yys(ii,jjj+1,kk+1,num1),yys(ii,jjj,kk+1,num1),yys(ii+1,jjj,kk,num1),yys(ii+1,jjj+1,kk,num1),yys(ii+1,jjj+1,kk+1,num1),yys(ii+1,jjj,kk+1,num1),&
                                    zzs(ii,jjj,kk,num1),zzs(ii,jjj+1,kk,num1),zzs(ii,jjj+1,kk+1,num1),zzs(ii,jjj,kk+1,num1),zzs(ii+1,jjj,kk,num1),zzs(ii+1,jjj+1,kk,num1),zzs(ii+1,jjj+1,kk+1,num1),zzs(ii+1,jjj,kk+1,num1),&
                                    x(i,j,k),x(i,j+1,k),x(i,j+1,k+1),x(i,j,k+1),x(i+1,j,k),x(i+1,j+1,k),x(i+1,j+1,k+1),x(i+1,j,k+1),&
                                    y(i,j,k),y(i,j+1,k),y(i,j+1,k+1),y(i,j,k+1),y(i+1,j,k),y(i+1,j+1,k),y(i+1,j+1,k+1),y(i+1,j,k+1),&
                                    z(i,j,k),z(i,j+1,k),z(i,j+1,k+1),z(i,j,k+1),z(i+1,j,k),z(i+1,j+1,k),z(i+1,j+1,k+1),z(i+1,j,k+1),overvv)
                                    v(j,k,ii,jjj,kk,l)=overvv/vv(j,k)
                            end if
                            if(yj/=j)then
                                jjj=yj
                                overvv=0.
                                call clipv(xxs(ii,jjj,kk,num1),xxs(ii,jjj+1,kk,num1),xxs(ii,jjj+1,kk+1,num1),xxs(ii,jjj,kk+1,num1),xxs(ii+1,jjj,kk,num1),xxs(ii+1,jjj+1,kk,num1),xxs(ii+1,jjj+1,kk+1,num1),xxs(ii+1,jjj,kk+1,num1),&
                                    yys(ii,jjj,kk,num1),yys(ii,jjj+1,kk,num1),yys(ii,jjj+1,kk+1,num1),yys(ii,jjj,kk+1,num1),yys(ii+1,jjj,kk,num1),yys(ii+1,jjj+1,kk,num1),yys(ii+1,jjj+1,kk+1,num1),yys(ii+1,jjj,kk+1,num1),&
                                    zzs(ii,jjj,kk,num1),zzs(ii,jjj+1,kk,num1),zzs(ii,jjj+1,kk+1,num1),zzs(ii,jjj,kk+1,num1),zzs(ii+1,jjj,kk,num1),zzs(ii+1,jjj+1,kk,num1),zzs(ii+1,jjj+1,kk+1,num1),zzs(ii+1,jjj,kk+1,num1),&
                                    x(i,j,k),x(i,j+1,k),x(i,j+1,k+1),x(i,j,k+1),x(i+1,j,k),x(i+1,j+1,k),x(i+1,j+1,k+1),x(i+1,j,k+1),&
                                    y(i,j,k),y(i,j+1,k),y(i,j+1,k+1),y(i,j,k+1),y(i+1,j,k),y(i+1,j+1,k),y(i+1,j+1,k+1),y(i+1,j,k+1),&
                                    z(i,j,k),z(i,j+1,k),z(i,j+1,k+1),z(i,j,k+1),z(i+1,j,k),z(i+1,j+1,k),z(i+1,j+1,k+1),z(i+1,j,k+1),overvv)
                                    v(j,k,ii,jjj,kk,l)=overvv/vv(j,k)
                            end if
                        end if
                        if(lll/=0 .and. overv==0.)then!决定了相邻叶排进程被裁剪循环的终止
                            exit
                        end if
                    end do
                    if(lll==0)then
                        kkk(ii)=1    !若当前进程某k与相邻叶排进程的所有k均不重叠，则退出，决定了当前进程裁剪循环的终止。2018.8.10新加的
                    end if
                end do
                if(iii==1 .and. kkk(1)==1 .and. kkk(2)==1)then
                    exit
                end if
            end do
        end do
    end do
    !各体积比求和不为1，进行后处理,按比例增加s(j,k,jj,kk,l)*(1.0-as)/as或扣除s(j,k,jj,kk,l)*(as-1.0)/as
    do j=1,ny
        do k=1,nz
            av=0.
            do l=1,ssum
                do jj=1,ny
                    do kk=1,nz
                        do ii=1,2
                            if(v(j,k,ii,jj,kk,l)/=0.)then
                                av=av+v(j,k,ii,jj,kk,l)
                            end if
                        end do
                    end do
                end do
            end do
            do l=1,ssum
                do jj=1,ny
                    do kk=1,nz
                        do ii=1,2
                            if(v(j,k,ii,jj,kk,l)/=0.)then
                                v(j,k,ii,jj,kk,l)=v(j,k,ii,jj,kk,l)/av
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do
    end subroutine overlapgrid
    
    subroutine overlapadj !确定重叠部分占<相邻进程>的jk范围，传递给相邻进程即bc
    use global   !与（myidl）的重叠部分对应《相邻进程》的小k，相邻进程（myidm）的中k，相邻进程(myidr)的大k。
    implicit none
    integer :: begink,endk,num
    
    bk=nz   !赋初值
    ek=1
    bkk=0
    ekk=0
    do l=1,ssum
        do j=1,ny
            do ii=1,2
                begink=0
                endk=0
                do kk=1,nz
                    jj=0
                    do k=1,nz
                        if(v(j,k,ii,j,kk,l)/=0.)then
                            jj=1
                            begink=kk
                            exit
                        end if
                    end do
                    if(jj==1)then
                        exit
                    end if
                end do
                do kk=nz,1,-1
                    jj=0
                    do k=1,nz
                        if(v(j,k,ii,j,kk,l)/=0.)then
                            jj=1
                            endk=kk
                            exit
                        end if
                    end do
                    if(jj==1)then
                        exit
                    end if
                end do
                if(begink<=bk(l) .and. begink/=0)then
                    bk(l)=begink
                end if
                if(endk>=ek(l) .and. endk/=0)then
                    ek(l)=endk
                end if
                if(bk(l)==1 .and. ek(l)==nz)then
                    exit
                end if
            end do
        end do
    end do
    do l=1,ssum   !bk，ek即重叠部分占相邻进程范围，从相邻进程接收数据后存放在该位置
        num=cli(l)!bkk，ekk即重叠部分占当前进程范围，当前进程按此范围发送数据
        call MPI_SENDRECV(bk(l),1,MPI_INTEGER,num,126,bkk(l),1,MPI_INTEGER,num,126,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(ek(l),1,MPI_INTEGER,num,127,ekk(l),1,MPI_INTEGER,num,127,MPI_COMM_WORLD,status,ierr)
    end do
    end subroutine overlapadj
   
    subroutine clip(x11,x22,x33,x44,x55,x66,x77,x88,y11,y22,y33,y44,y55,y66,y77,y88,z11,z22,z33,z44,z55,z66,z77,z88,&
                    x011,x022,x033,x044,x055,x066,x077,x088,y011,y022,y033,y044,y055,y066,y077,y088,z011,z022,z033,z044,z055,z066,z077,z088,i4)
    implicit none
    real(8),intent(in) :: x11,x22,x33,x44,x55,x66,x77,x88,y11,y22,y33,y44,y55,y66,y77,y88,z11,z22,z33,z44,z55,z66,z77,z88
    real(8),intent(in) :: x011,x022,x033,x044,x055,x066,x077,x088,y011,y022,y033,y044,y055,y066,y077,y088,z011,z022,z033,z044,z055,z066,z077,z088
    integer,intent(out) ::i4
    real(8) :: bjx(0:1,0:1,0:1),bjy(0:1,0:1,0:1),bjz(0:1,0:1,0:1),jdx(0:1,0:1,0:1),jdy(0:1,0:1,0:1),jdz(0:1,0:1,0:1),sbjx(6,4),sbjy(6,4),sbjz(6,4),sjdx(6,4),sjdy(6,4),sjdz(6,4)
    real(8) :: shax(6),shay(6),shaz(6),sjdxn(6),sjdyn(6),sjdzn(6),osx(20,20),osy(20,20),osz(20,20),x(50,50),y(50,50),z(50,50),xj(50),yj(50),zj(50),xjj(20),yjj(20),zjj(20),osx0(20),osy0(20),osz0(20)
    integer :: i,j,k,ii,os,oss,sum,num,nnum,j1,j2,kk,sum0,l
    integer :: sa(20,20,6),sn(20),od(20),odd(20)
    real(8) :: jdx0,jdy0,jdz0,x0,y0,z0,xx,yy,zz,rr,x21,y21,z21,x23,y23,z23,rrx,rry,rrz,dd
    !剪刀体的8个顶点,6个面(4个顶点)
    jdx(0,0,0)=x11
    jdx(0,1,0)=x22
    jdx(0,1,1)=x33
    jdx(0,0,1)=x44
    jdx(1,0,0)=x55
    jdx(1,1,0)=x66
    jdx(1,1,1)=x77
    jdx(1,0,1)=x88
    jdy(0,0,0)=y11
    jdy(0,1,0)=y22
    jdy(0,1,1)=y33
    jdy(0,0,1)=y44
    jdy(1,0,0)=y55
    jdy(1,1,0)=y66
    jdy(1,1,1)=y77
    jdy(1,0,1)=y88
    jdz(0,0,0)=z11
    jdz(0,1,0)=z22
    jdz(0,1,1)=z33
    jdz(0,0,1)=z44
    jdz(1,0,0)=z55
    jdz(1,1,0)=z66
    jdz(1,1,1)=z77
    jdz(1,0,1)=z88
    sjdx(1,1)=jdx(0,0,0)   !6个面,每个面上的4个顶点
    sjdy(1,1)=jdy(0,0,0)   !左面i=0,面1
    sjdz(1,1)=jdz(0,0,0)
    sjdx(1,2)=jdx(0,1,0)
    sjdy(1,2)=jdy(0,1,0)
    sjdz(1,2)=jdz(0,1,0)
    sjdx(1,3)=jdx(0,1,1)
    sjdy(1,3)=jdy(0,1,1)
    sjdz(1,3)=jdz(0,1,1)
    sjdx(1,4)=jdx(0,0,1)
    sjdy(1,4)=jdy(0,0,1)
    sjdz(1,4)=jdz(0,0,1)
    sjdx(2,1)=jdx(1,0,0)   !右面i=1,面2
    sjdy(2,1)=jdy(1,0,0)
    sjdz(2,1)=jdz(1,0,0)
    sjdx(2,2)=jdx(1,1,0)
    sjdy(2,2)=jdy(1,1,0)
    sjdz(2,2)=jdz(1,1,0)
    sjdx(2,3)=jdx(1,1,1)
    sjdy(2,3)=jdy(1,1,1)
    sjdz(2,3)=jdz(1,1,1)
    sjdx(2,4)=jdx(1,0,1)
    sjdy(2,4)=jdy(1,0,1)
    sjdz(2,4)=jdz(1,0,1)
    sjdx(3,1)=jdx(0,0,0)   !下面j=0,面3
    sjdy(3,1)=jdy(0,0,0)
    sjdz(3,1)=jdz(0,0,0)
    sjdx(3,2)=jdx(1,0,0)
    sjdy(3,2)=jdy(1,0,0)
    sjdz(3,2)=jdz(1,0,0)
    sjdx(3,3)=jdx(1,0,1)
    sjdy(3,3)=jdy(1,0,1)
    sjdz(3,3)=jdz(1,0,1)
    sjdx(3,4)=jdx(0,0,1)
    sjdy(3,4)=jdy(0,0,1)
    sjdz(3,4)=jdz(0,0,1)
    sjdx(4,1)=jdx(0,1,0)   !上面j=1,面4
    sjdy(4,1)=jdy(0,1,0)
    sjdz(4,1)=jdz(0,1,0)
    sjdx(4,2)=jdx(1,1,0)
    sjdy(4,2)=jdy(1,1,0)
    sjdz(4,2)=jdz(1,1,0)
    sjdx(4,3)=jdx(1,1,1)
    sjdy(4,3)=jdy(1,1,1)
    sjdz(4,3)=jdz(1,1,1)
    sjdx(4,4)=jdx(0,1,1)
    sjdy(4,4)=jdy(0,1,1)
    sjdz(4,4)=jdz(0,1,1)
    sjdx(5,1)=jdx(0,0,0)   !后面k=0,面5
    sjdy(5,1)=jdy(0,0,0)
    sjdz(5,1)=jdz(0,0,0)
    sjdx(5,2)=jdx(1,0,0)
    sjdy(5,2)=jdy(1,0,0)
    sjdz(5,2)=jdz(1,0,0)
    sjdx(5,3)=jdx(1,1,0)
    sjdy(5,3)=jdy(1,1,0)
    sjdz(5,3)=jdz(1,1,0)
    sjdx(5,4)=jdx(0,1,0)
    sjdy(5,4)=jdy(0,1,0)
    sjdz(5,4)=jdz(0,1,0)
    sjdx(6,1)=jdx(0,0,1)   !前面k=1,面6
    sjdy(6,1)=jdy(0,0,1)
    sjdz(6,1)=jdz(0,0,1)
    sjdx(6,2)=jdx(1,0,1)
    sjdy(6,2)=jdy(1,0,1)
    sjdz(6,2)=jdz(1,0,1)
    sjdx(6,3)=jdx(1,1,1)
    sjdy(6,3)=jdy(1,1,1)
    sjdz(6,3)=jdz(1,1,1)
    sjdx(6,4)=jdx(0,1,1)
    sjdy(6,4)=jdy(0,1,1)
    sjdz(6,4)=jdz(0,1,1)
    !被剪体的8个顶点,6个面(4个顶点)
    bjx(0,0,0)=x011   !8个顶点
    bjx(0,1,0)=x022
    bjx(0,1,1)=x033
    bjx(0,0,1)=x044
    bjx(1,0,0)=x055
    bjx(1,1,0)=x066
    bjx(1,1,1)=x077
    bjx(1,0,1)=x088
    bjy(0,0,0)=y011
    bjy(0,1,0)=y022
    bjy(0,1,1)=y033
    bjy(0,0,1)=y044
    bjy(1,0,0)=y055
    bjy(1,1,0)=y066
    bjy(1,1,1)=y077
    bjy(1,0,1)=y088
    bjz(0,0,0)=z011
    bjz(0,1,0)=z022
    bjz(0,1,1)=z033
    bjz(0,0,1)=z044
    bjz(1,0,0)=z055
    bjz(1,1,0)=z066
    bjz(1,1,1)=z077
    bjz(1,0,1)=z088
    sbjx(1,1)=bjx(0,0,0)   !6个面,每个面上的4个顶点
    sbjy(1,1)=bjy(0,0,0)   !左面i=0,面1
    sbjz(1,1)=bjz(0,0,0)
    sbjx(1,2)=bjx(0,1,0)
    sbjy(1,2)=bjy(0,1,0)
    sbjz(1,2)=bjz(0,1,0)
    sbjx(1,3)=bjx(0,1,1)
    sbjy(1,3)=bjy(0,1,1)
    sbjz(1,3)=bjz(0,1,1)
    sbjx(1,4)=bjx(0,0,1)
    sbjy(1,4)=bjy(0,0,1)
    sbjz(1,4)=bjz(0,0,1)
    sbjx(2,1)=bjx(1,0,0)   !右面i=1,面2
    sbjy(2,1)=bjy(1,0,0)
    sbjz(2,1)=bjz(1,0,0)
    sbjx(2,2)=bjx(1,1,0)
    sbjy(2,2)=bjy(1,1,0)
    sbjz(2,2)=bjz(1,1,0)
    sbjx(2,3)=bjx(1,1,1)
    sbjy(2,3)=bjy(1,1,1)
    sbjz(2,3)=bjz(1,1,1)
    sbjx(2,4)=bjx(1,0,1)
    sbjy(2,4)=bjy(1,0,1)
    sbjz(2,4)=bjz(1,0,1)
    sbjx(3,1)=bjx(0,0,0)   !下面j=0,面3
    sbjy(3,1)=bjy(0,0,0)
    sbjz(3,1)=bjz(0,0,0)
    sbjx(3,2)=bjx(1,0,0)
    sbjy(3,2)=bjy(1,0,0)
    sbjz(3,2)=bjz(1,0,0)
    sbjx(3,3)=bjx(1,0,1)
    sbjy(3,3)=bjy(1,0,1)
    sbjz(3,3)=bjz(1,0,1)
    sbjx(3,4)=bjx(0,0,1)
    sbjy(3,4)=bjy(0,0,1)
    sbjz(3,4)=bjz(0,0,1)
    sbjx(4,1)=bjx(0,1,0)   !上面j=1,面4
    sbjy(4,1)=bjy(0,1,0)
    sbjz(4,1)=bjz(0,1,0)
    sbjx(4,2)=bjx(1,1,0)
    sbjy(4,2)=bjy(1,1,0)
    sbjz(4,2)=bjz(1,1,0)
    sbjx(4,3)=bjx(1,1,1)
    sbjy(4,3)=bjy(1,1,1)
    sbjz(4,3)=bjz(1,1,1)
    sbjx(4,4)=bjx(0,1,1)
    sbjy(4,4)=bjy(0,1,1)
    sbjz(4,4)=bjz(0,1,1)
    sbjx(5,1)=bjx(0,0,0)   !后面k=0,面5
    sbjy(5,1)=bjy(0,0,0)
    sbjz(5,1)=bjz(0,0,0)
    sbjx(5,2)=bjx(1,0,0)
    sbjy(5,2)=bjy(1,0,0)
    sbjz(5,2)=bjz(1,0,0)
    sbjx(5,3)=bjx(1,1,0)
    sbjy(5,3)=bjy(1,1,0)
    sbjz(5,3)=bjz(1,1,0)
    sbjx(5,4)=bjx(0,1,0)
    sbjy(5,4)=bjy(0,1,0)
    sbjz(5,4)=bjz(0,1,0)
    sbjx(6,1)=bjx(0,0,1)   !前面k=1,面6
    sbjy(6,1)=bjy(0,0,1)
    sbjz(6,1)=bjz(0,0,1)
    sbjx(6,2)=bjx(1,0,1)
    sbjy(6,2)=bjy(1,0,1)
    sbjz(6,2)=bjz(1,0,1)
    sbjx(6,3)=bjx(1,1,1)
    sbjy(6,3)=bjy(1,1,1)
    sbjz(6,3)=bjz(1,1,1)
    sbjx(6,4)=bjx(0,1,1)
    sbjy(6,4)=bjy(0,1,1)
    sbjz(6,4)=bjz(0,1,1)
    jdx0=0.125d0*(jdx(0,0,0)+jdx(0,1,0)+jdx(0,1,1)+jdx(0,0,1)+jdx(1,0,0)+jdx(1,1,0)+jdx(1,1,1)+jdx(1,0,1))!剪刀体中心坐标
    jdy0=0.125d0*(jdy(0,0,0)+jdy(0,1,0)+jdy(0,1,1)+jdy(0,0,1)+jdy(1,0,0)+jdy(1,1,0)+jdy(1,1,1)+jdy(1,0,1))
    jdz0=0.125d0*(jdz(0,0,0)+jdz(0,1,0)+jdz(0,1,1)+jdz(0,0,1)+jdz(1,0,0)+jdz(1,1,0)+jdz(1,1,1)+jdz(1,0,1))
    do i=1,6
        call shadow(sjdx(i,1),sjdy(i,1),sjdz(i,1),sjdx(i,2),sjdy(i,2),sjdz(i,2),sjdx(i,3),sjdy(i,3),sjdz(i,3),jdx0,jdy0,jdz0,shax(i),shay(i),shaz(i))!体中心到面的投影点坐标
        x21=sjdx(i,1)-sjdx(i,2)
        y21=sjdy(i,1)-sjdy(i,2)
        z21=sjdz(i,1)-sjdz(i,2)
        x23=sjdx(i,3)-sjdx(i,2)
        y23=sjdy(i,3)-sjdy(i,2)
        z23=sjdz(i,3)-sjdz(i,2)
        sjdxn(i)=y21*z23-z21*y23  !剪刀体每个面的法向量
        sjdyn(i)=z21*x23-x21*z23
        sjdzn(i)=x21*y23-y21*x23
    end do
    do ii=1,6!剪刀体6个面所在的无界面进行裁剪,外循环
        oss=0
        odd=0
        x=0.
        y=0.
        z=0.
        if(ii==1)then
            oss=6
            do i=1,oss
                odd(i)=4
                do j=1,odd(i)
                    x(i,j)=sbjx(i,j)
                    y(i,j)=sbjy(i,j)
                    z(i,j)=sbjz(i,j)
                end do
            end do
        else
            if(os==0)then
                exit
            else
                oss=os
                do i=1,oss
                    odd(i)=od(i)
                    do j=1,odd(i)
                        x(i,j)=osx(i,j)
                        y(i,j)=osy(i,j)
                        z(i,j)=osz(i,j)
                    end do
                end do
            end if
        end if
        sa=3
        sn=3
        do i=1,oss
            !判断被剪体当前所有面上所有顶点在剪刀体当前面的内侧0\面上2\外侧1
            do j=1,odd(i)
                rr=(shax(ii)-x(i,j))*(shax(ii)-jdx0)+(shay(ii)-y(i,j))*(shay(ii)-jdy0)+(shaz(ii)-z(i,j))*(shaz(ii)-jdz0)
                if(rr>0.)then!夹角小于90度,被剪顶点与剪刀重心一样在该面内
                    sa(i,j,ii)=0
                else if(rr==0.)then!夹角等于90度,被剪顶点刚好在面上
                    sa(i,j,ii)=2
                else!夹角大于90度,被剪顶点与剪刀重心在该面两侧
                    sa(i,j,ii)=1
                end if
            end do
            !判断被剪体当前所有面与剪刀体当前面是否平行
            x21=x(i,1)-x(i,2)
            y21=y(i,1)-y(i,2)
            z21=z(i,1)-z(i,2)
            x23=x(i,3)-x(i,2)
            y23=y(i,3)-y(i,2)
            z23=z(i,3)-z(i,2)
            x0=y21*z23-z21*y23  !被剪体每个面的法向量
            y0=z21*x23-x21*z23
            z0=x21*y23-y21*x23
            xx=sjdyn(ii)*z0-sjdzn(ii)*y0!剪刀体与被剪体两两面法向量叉乘
            yy=sjdzn(ii)*x0-sjdxn(ii)*z0
            zz=sjdxn(ii)*y0-sjdyn(ii)*x0
            rr=xx*xx+yy*yy+zz*zz
            if(rr==0.)then
                sn(i)=0 !两面平行
            else
                sn(i)=1 !两面不平行
            end if
        end do
        
        os=0   !裁剪后得到的面个数
        od=0   !裁剪后每个面上的点个数
        num=0  !所有面上裁剪得到的交点总数
        xj=0.
        yj=0.
        zj=0.
        do i=1,oss!被剪体的有界面,内循环
            if(sn(i)==0)then!面平行，判断其中一个点是面内上外侧即可
                if(sa(i,1,ii)==1)then  !在面外侧,全不要
                else    !在面上、内侧,全要
                    os=os+1
                    od(os)=odd(i)
                    do j=1,od(os)
                        osx(os,j)=x(i,j)
                        osy(os,j)=y(i,j)
                        osz(os,j)=z(i,j)
                    end do
                end if
            else  !两面不平行，判断面所有顶点
                k=0
                do j=1,odd(i)
                    if(sa(i,j,ii)/=0)then!所有顶点都在面外或面上，都不要
                        k=k+1
                    end if
                end do
                if(k==odd(i))then
                else
                    os=os+1
                    sum=0
                    do j=1,odd(i)
                        j1=j
                        j2=j+1
                        if(j==odd(i))then
                            j2=1
                        end if
                        if(sa(i,j1,ii)==0 .and. sa(i,j2,ii)==0)then
                            osx(os,sum+1)=x(i,j1)
                            osy(os,sum+1)=y(i,j1)
                            osz(os,sum+1)=z(i,j1)
                            osx(os,sum+2)=x(i,j2)
                            osy(os,sum+2)=y(i,j2)
                            osz(os,sum+2)=z(i,j2)
                            sum=sum+2
                        else if(sa(i,j1,ii)==0 .and. sa(i,j2,ii)==1)then
                            call interls(sjdx(ii,1),sjdy(ii,1),sjdz(ii,1),sjdx(ii,2),sjdy(ii,2),sjdz(ii,2),sjdx(ii,3),sjdy(ii,3),sjdz(ii,3),x(i,j1),y(i,j1),z(i,j1),x(i,j2),y(i,j2),z(i,j2),xx,yy,zz)
                            osx(os,sum+1)=x(i,j1)
                            osy(os,sum+1)=y(i,j1)
                            osz(os,sum+1)=z(i,j1)
                            osx(os,sum+2)=xx
                            osy(os,sum+2)=yy
                            osz(os,sum+2)=zz
                            sum=sum+2
                            xj(num+1)=xx
                            yj(num+1)=yy
                            zj(num+1)=zz
                            num=num+1
                        else if(sa(i,j1,ii)==0 .and. sa(i,j2,ii)==2)then
                            osx(os,sum+1)=x(i,j1)
                            osy(os,sum+1)=y(i,j1)
                            osz(os,sum+1)=z(i,j1)
                            osx(os,sum+2)=x(i,j2)
                            osy(os,sum+2)=y(i,j2)
                            osz(os,sum+2)=z(i,j2)
                            sum=sum+2
                            xj(num+1)=x(i,j2)
                            yj(num+1)=y(i,j2)
                            zj(num+1)=z(i,j2)
                            num=num+1
                        else if(sa(i,j1,ii)==1 .and. sa(i,j2,ii)==0)then
                            call interls(sjdx(ii,1),sjdy(ii,1),sjdz(ii,1),sjdx(ii,2),sjdy(ii,2),sjdz(ii,2),sjdx(ii,3),sjdy(ii,3),sjdz(ii,3),x(i,j1),y(i,j1),z(i,j1),x(i,j2),y(i,j2),z(i,j2),xx,yy,zz)
                            osx(os,sum+1)=xx
                            osy(os,sum+1)=yy
                            osz(os,sum+1)=zz
                            osx(os,sum+2)=x(i,j2)
                            osy(os,sum+2)=y(i,j2)
                            osz(os,sum+2)=z(i,j2)
                            sum=sum+2
                            xj(num+1)=xx
                            yj(num+1)=yy
                            zj(num+1)=zz
                            num=num+1
                        else if(sa(i,j1,ii)==1 .and. sa(i,j2,ii)==1)then
                        else if(sa(i,j1,ii)==1 .and. sa(i,j2,ii)==2)then
                        else if(sa(i,j1,ii)==2 .and. sa(i,j2,ii)==0)then
                            osx(os,sum+1)=x(i,j1)
                            osy(os,sum+1)=y(i,j1)
                            osz(os,sum+1)=z(i,j1)
                            osx(os,sum+2)=x(i,j2)
                            osy(os,sum+2)=y(i,j2)
                            osz(os,sum+2)=z(i,j2)
                            sum=sum+2
                            xj(num+1)=x(i,j1)
                            yj(num+1)=y(i,j1)
                            zj(num+1)=z(i,j1)
                            num=num+1
                        else if(sa(i,j1,ii)==2 .and. sa(i,j2,ii)==1)then
                        else if(sa(i,j1,ii)==2 .and. sa(i,j2,ii)==2)then
                            osx(os,sum+1)=x(i,j1)
                            osy(os,sum+1)=y(i,j1)
                            osz(os,sum+1)=z(i,j1)
                            osx(os,sum+2)=x(i,j2)
                            osy(os,sum+2)=y(i,j2)
                            osz(os,sum+2)=z(i,j2)
                            sum=sum+2
                            xj(num+1)=x(i,j1)
                            yj(num+1)=y(i,j1)
                            zj(num+1)=z(i,j1)
                            xj(num+2)=x(i,j2)
                            yj(num+2)=y(i,j2)
                            zj(num+2)=z(i,j2)
                            num=num+2
                        end if
                    end do
                    if(sum>1)then
                        sum0=1
                        osx0(1)=osx(os,1)
                        osy0(1)=osy(os,1)
                        osz0(1)=osz(os,1)
                        do k=2,sum
                            l=0
                            do kk=1,sum0
                                if(abs(osx0(kk)-osx(os,k))<1.d-12 .and. abs(osy0(kk)-osy(os,k))<1.d-12 .and. abs(osz0(kk)-osz(os,k))<1.d-12)then
                                    l=1
                                    exit
                                end if
                            end do
                            if(l==0)then !该点与记录里的所有点都不重合，则记录该点
                                sum0=sum0+1
                                osx0(sum0)=osx(os,k)
                                osy0(sum0)=osy(os,k)
                                osz0(sum0)=osz(os,k)
                            end if
                        end do
                        od(os)=sum0
                        osx(os,1:sum0)=osx0(1:sum0)
                        osy(os,1:sum0)=osy0(1:sum0)
                        osz(os,1:sum0)=osz0(1:sum0)
                    end if
                end if
            end if
        end do
        if(num>0)then!因裁剪出来的新界面
            i=0
            call newsec(num,xj(1:num),yj(1:num),zj(1:num),i,nnum,xjj(1:nnum),yjj(1:nnum),zjj(1:nnum))
            if(i==0)then
            else!得到新面，还要判断与某os面是否共面，一旦有就不要新面
                sum=0
                x21=xjj(1)-xjj(2)
                y21=yjj(1)-yjj(2)
                z21=zjj(1)-zjj(2)
                x23=xjj(3)-xjj(2)
                y23=yjj(3)-yjj(2)
                z23=zjj(3)-zjj(2)
                x0=y21*z23-z21*y23
                y0=z21*x23-x21*z23
                z0=x21*y23-y21*x23
                do i=1,os
                    x21=osx(i,1)-osx(i,2)
                    y21=osy(i,1)-osy(i,2)
                    z21=osz(i,1)-osz(i,2)
                    x23=osx(i,3)-osx(i,2)
                    y23=osy(i,3)-osy(i,2)
                    z23=osz(i,3)-osz(i,2)
                    xx=y21*z23-z21*y23
                    yy=z21*x23-x21*z23
                    zz=x21*y23-y21*x23
                    rrx=yy*z0-zz*y0
                    rry=zz*x0-xx*z0
                    rrz=xx*y0-yy*x0
                    rr=rrx**2+rry**2+rrz**2
                    if(rr==0.)then!两面平行
                        call shadow(xjj(1),yjj(1),zjj(1),xjj(2),yjj(2),zjj(2),xjj(3),yjj(3),zjj(3),osx(i,1),osy(i,1),osz(i,1),xx,yy,zz)
                        dd=(osx(i,1)-xx)**2+(osy(i,1)-yy)**2+(osz(i,1)-zz)**2
                        if(dd==0.)then!两面重合
                            sum=1
                            exit
                        end if
                    end if
                end do
                if(sum==0)then
                    os=os+1
                    od(os)=nnum
                    do j=1,od(os)
                        osx(os,j)=xjj(j)
                        osy(os,j)=yjj(j)
                        osz(os,j)=zjj(j)
                    end do
                end if
            end if
        end if
    end do
    i4=0
    if(os>=4)then
        i4=1
    end if
    end subroutine clip
    
    subroutine clipv(x11,x22,x33,x44,x55,x66,x77,x88,y11,y22,y33,y44,y55,y66,y77,y88,z11,z22,z33,z44,z55,z66,z77,z88,&
                    x011,x022,x033,x044,x055,x066,x077,x088,y011,y022,y033,y044,y055,y066,y077,y088,z011,z022,z033,z044,z055,z066,z077,z088,v)
    implicit none
    real(8),intent(in) :: x11,x22,x33,x44,x55,x66,x77,x88,y11,y22,y33,y44,y55,y66,y77,y88,z11,z22,z33,z44,z55,z66,z77,z88
    real(8),intent(in) :: x011,x022,x033,x044,x055,x066,x077,x088,y011,y022,y033,y044,y055,y066,y077,y088,z011,z022,z033,z044,z055,z066,z077,z088
    real(8),intent(out) ::v
    real(8) :: bjx(0:1,0:1,0:1),bjy(0:1,0:1,0:1),bjz(0:1,0:1,0:1),jdx(0:1,0:1,0:1),jdy(0:1,0:1,0:1),jdz(0:1,0:1,0:1),sbjx(6,4),sbjy(6,4),sbjz(6,4),sjdx(6,4),sjdy(6,4),sjdz(6,4)
    real(8) :: shax(6),shay(6),shaz(6),sjdxn(6),sjdyn(6),sjdzn(6),osx(20,20),osy(20,20),osz(20,20),x(50,50),y(50,50),z(50,50),xj(50),yj(50),zj(50),xjj(20),yjj(20),zjj(20)
    real(8) :: zdx(20),zdy(20),zdz(20),osx0(20),osy0(20),osz0(20)
    integer :: i,j,k,l,ii,os,oss,sum,num,nnum,j1,j2,kk,sum0
    integer :: sa(20,20,6),sn(20),od(20),odd(20)
    real(8) :: jdx0,jdy0,jdz0,x0,y0,z0,xx,yy,zz,rr,x21,y21,z21,x23,y23,z23,rrx,rry,rrz,dd,rx,ry,rz,hh
    !剪刀体的8个顶点,6个面(4个顶点)
    jdx(0,0,0)=x11
    jdx(0,1,0)=x22
    jdx(0,1,1)=x33
    jdx(0,0,1)=x44
    jdx(1,0,0)=x55
    jdx(1,1,0)=x66
    jdx(1,1,1)=x77
    jdx(1,0,1)=x88
    jdy(0,0,0)=y11
    jdy(0,1,0)=y22
    jdy(0,1,1)=y33
    jdy(0,0,1)=y44
    jdy(1,0,0)=y55
    jdy(1,1,0)=y66
    jdy(1,1,1)=y77
    jdy(1,0,1)=y88
    jdz(0,0,0)=z11
    jdz(0,1,0)=z22
    jdz(0,1,1)=z33
    jdz(0,0,1)=z44
    jdz(1,0,0)=z55
    jdz(1,1,0)=z66
    jdz(1,1,1)=z77
    jdz(1,0,1)=z88
    sjdx(1,1)=jdx(0,0,0)   !6个面,每个面上的4个顶点
    sjdy(1,1)=jdy(0,0,0)   !左面i=0,面1
    sjdz(1,1)=jdz(0,0,0)
    sjdx(1,2)=jdx(0,1,0)
    sjdy(1,2)=jdy(0,1,0)
    sjdz(1,2)=jdz(0,1,0)
    sjdx(1,3)=jdx(0,1,1)
    sjdy(1,3)=jdy(0,1,1)
    sjdz(1,3)=jdz(0,1,1)
    sjdx(1,4)=jdx(0,0,1)
    sjdy(1,4)=jdy(0,0,1)
    sjdz(1,4)=jdz(0,0,1)
    sjdx(2,1)=jdx(1,0,0)   !右面i=1,面2
    sjdy(2,1)=jdy(1,0,0)
    sjdz(2,1)=jdz(1,0,0)
    sjdx(2,2)=jdx(1,1,0)
    sjdy(2,2)=jdy(1,1,0)
    sjdz(2,2)=jdz(1,1,0)
    sjdx(2,3)=jdx(1,1,1)
    sjdy(2,3)=jdy(1,1,1)
    sjdz(2,3)=jdz(1,1,1)
    sjdx(2,4)=jdx(1,0,1)
    sjdy(2,4)=jdy(1,0,1)
    sjdz(2,4)=jdz(1,0,1)
    sjdx(3,1)=jdx(0,0,0)   !下面j=0,面3
    sjdy(3,1)=jdy(0,0,0)
    sjdz(3,1)=jdz(0,0,0)
    sjdx(3,2)=jdx(1,0,0)
    sjdy(3,2)=jdy(1,0,0)
    sjdz(3,2)=jdz(1,0,0)
    sjdx(3,3)=jdx(1,0,1)
    sjdy(3,3)=jdy(1,0,1)
    sjdz(3,3)=jdz(1,0,1)
    sjdx(3,4)=jdx(0,0,1)
    sjdy(3,4)=jdy(0,0,1)
    sjdz(3,4)=jdz(0,0,1)
    sjdx(4,1)=jdx(0,1,0)   !上面j=1,面4
    sjdy(4,1)=jdy(0,1,0)
    sjdz(4,1)=jdz(0,1,0)
    sjdx(4,2)=jdx(1,1,0)
    sjdy(4,2)=jdy(1,1,0)
    sjdz(4,2)=jdz(1,1,0)
    sjdx(4,3)=jdx(1,1,1)
    sjdy(4,3)=jdy(1,1,1)
    sjdz(4,3)=jdz(1,1,1)
    sjdx(4,4)=jdx(0,1,1)
    sjdy(4,4)=jdy(0,1,1)
    sjdz(4,4)=jdz(0,1,1)
    sjdx(5,1)=jdx(0,0,0)   !后面k=0,面5
    sjdy(5,1)=jdy(0,0,0)
    sjdz(5,1)=jdz(0,0,0)
    sjdx(5,2)=jdx(1,0,0)
    sjdy(5,2)=jdy(1,0,0)
    sjdz(5,2)=jdz(1,0,0)
    sjdx(5,3)=jdx(1,1,0)
    sjdy(5,3)=jdy(1,1,0)
    sjdz(5,3)=jdz(1,1,0)
    sjdx(5,4)=jdx(0,1,0)
    sjdy(5,4)=jdy(0,1,0)
    sjdz(5,4)=jdz(0,1,0)
    sjdx(6,1)=jdx(0,0,1)   !前面k=1,面6
    sjdy(6,1)=jdy(0,0,1)
    sjdz(6,1)=jdz(0,0,1)
    sjdx(6,2)=jdx(1,0,1)
    sjdy(6,2)=jdy(1,0,1)
    sjdz(6,2)=jdz(1,0,1)
    sjdx(6,3)=jdx(1,1,1)
    sjdy(6,3)=jdy(1,1,1)
    sjdz(6,3)=jdz(1,1,1)
    sjdx(6,4)=jdx(0,1,1)
    sjdy(6,4)=jdy(0,1,1)
    sjdz(6,4)=jdz(0,1,1)
    !被剪体的8个顶点,6个面(4个顶点)
    bjx(0,0,0)=x011   !8个顶点
    bjx(0,1,0)=x022
    bjx(0,1,1)=x033
    bjx(0,0,1)=x044
    bjx(1,0,0)=x055
    bjx(1,1,0)=x066
    bjx(1,1,1)=x077
    bjx(1,0,1)=x088
    bjy(0,0,0)=y011
    bjy(0,1,0)=y022
    bjy(0,1,1)=y033
    bjy(0,0,1)=y044
    bjy(1,0,0)=y055
    bjy(1,1,0)=y066
    bjy(1,1,1)=y077
    bjy(1,0,1)=y088
    bjz(0,0,0)=z011
    bjz(0,1,0)=z022
    bjz(0,1,1)=z033
    bjz(0,0,1)=z044
    bjz(1,0,0)=z055
    bjz(1,1,0)=z066
    bjz(1,1,1)=z077
    bjz(1,0,1)=z088
    sbjx(1,1)=bjx(0,0,0)   !6个面,每个面上的4个顶点
    sbjy(1,1)=bjy(0,0,0)   !左面i=0,面1
    sbjz(1,1)=bjz(0,0,0)
    sbjx(1,2)=bjx(0,1,0)
    sbjy(1,2)=bjy(0,1,0)
    sbjz(1,2)=bjz(0,1,0)
    sbjx(1,3)=bjx(0,1,1)
    sbjy(1,3)=bjy(0,1,1)
    sbjz(1,3)=bjz(0,1,1)
    sbjx(1,4)=bjx(0,0,1)
    sbjy(1,4)=bjy(0,0,1)
    sbjz(1,4)=bjz(0,0,1)
    sbjx(2,1)=bjx(1,0,0)   !右面i=1,面2
    sbjy(2,1)=bjy(1,0,0)
    sbjz(2,1)=bjz(1,0,0)
    sbjx(2,2)=bjx(1,1,0)
    sbjy(2,2)=bjy(1,1,0)
    sbjz(2,2)=bjz(1,1,0)
    sbjx(2,3)=bjx(1,1,1)
    sbjy(2,3)=bjy(1,1,1)
    sbjz(2,3)=bjz(1,1,1)
    sbjx(2,4)=bjx(1,0,1)
    sbjy(2,4)=bjy(1,0,1)
    sbjz(2,4)=bjz(1,0,1)
    sbjx(3,1)=bjx(0,0,0)   !下面j=0,面3
    sbjy(3,1)=bjy(0,0,0)
    sbjz(3,1)=bjz(0,0,0)
    sbjx(3,2)=bjx(1,0,0)
    sbjy(3,2)=bjy(1,0,0)
    sbjz(3,2)=bjz(1,0,0)
    sbjx(3,3)=bjx(1,0,1)
    sbjy(3,3)=bjy(1,0,1)
    sbjz(3,3)=bjz(1,0,1)
    sbjx(3,4)=bjx(0,0,1)
    sbjy(3,4)=bjy(0,0,1)
    sbjz(3,4)=bjz(0,0,1)
    sbjx(4,1)=bjx(0,1,0)   !上面j=1,面4
    sbjy(4,1)=bjy(0,1,0)
    sbjz(4,1)=bjz(0,1,0)
    sbjx(4,2)=bjx(1,1,0)
    sbjy(4,2)=bjy(1,1,0)
    sbjz(4,2)=bjz(1,1,0)
    sbjx(4,3)=bjx(1,1,1)
    sbjy(4,3)=bjy(1,1,1)
    sbjz(4,3)=bjz(1,1,1)
    sbjx(4,4)=bjx(0,1,1)
    sbjy(4,4)=bjy(0,1,1)
    sbjz(4,4)=bjz(0,1,1)
    sbjx(5,1)=bjx(0,0,0)   !后面k=0,面5
    sbjy(5,1)=bjy(0,0,0)
    sbjz(5,1)=bjz(0,0,0)
    sbjx(5,2)=bjx(1,0,0)
    sbjy(5,2)=bjy(1,0,0)
    sbjz(5,2)=bjz(1,0,0)
    sbjx(5,3)=bjx(1,1,0)
    sbjy(5,3)=bjy(1,1,0)
    sbjz(5,3)=bjz(1,1,0)
    sbjx(5,4)=bjx(0,1,0)
    sbjy(5,4)=bjy(0,1,0)
    sbjz(5,4)=bjz(0,1,0)
    sbjx(6,1)=bjx(0,0,1)   !前面k=1,面6
    sbjy(6,1)=bjy(0,0,1)
    sbjz(6,1)=bjz(0,0,1)
    sbjx(6,2)=bjx(1,0,1)
    sbjy(6,2)=bjy(1,0,1)
    sbjz(6,2)=bjz(1,0,1)
    sbjx(6,3)=bjx(1,1,1)
    sbjy(6,3)=bjy(1,1,1)
    sbjz(6,3)=bjz(1,1,1)
    sbjx(6,4)=bjx(0,1,1)
    sbjy(6,4)=bjy(0,1,1)
    sbjz(6,4)=bjz(0,1,1)
    jdx0=0.125d0*(jdx(0,0,0)+jdx(0,1,0)+jdx(0,1,1)+jdx(0,0,1)+jdx(1,0,0)+jdx(1,1,0)+jdx(1,1,1)+jdx(1,0,1))!剪刀体中心坐标
    jdy0=0.125d0*(jdy(0,0,0)+jdy(0,1,0)+jdy(0,1,1)+jdy(0,0,1)+jdy(1,0,0)+jdy(1,1,0)+jdy(1,1,1)+jdy(1,0,1))
    jdz0=0.125d0*(jdz(0,0,0)+jdz(0,1,0)+jdz(0,1,1)+jdz(0,0,1)+jdz(1,0,0)+jdz(1,1,0)+jdz(1,1,1)+jdz(1,0,1))
    do i=1,6
        call shadow(sjdx(i,1),sjdy(i,1),sjdz(i,1),sjdx(i,2),sjdy(i,2),sjdz(i,2),sjdx(i,3),sjdy(i,3),sjdz(i,3),jdx0,jdy0,jdz0,shax(i),shay(i),shaz(i))!体中心到面的投影点坐标
        x21=sjdx(i,1)-sjdx(i,2)
        y21=sjdy(i,1)-sjdy(i,2)
        z21=sjdz(i,1)-sjdz(i,2)
        x23=sjdx(i,3)-sjdx(i,2)
        y23=sjdy(i,3)-sjdy(i,2)
        z23=sjdz(i,3)-sjdz(i,2)
        sjdxn(i)=y21*z23-z21*y23  !剪刀体每个面的法向量
        sjdyn(i)=z21*x23-x21*z23
        sjdzn(i)=x21*y23-y21*x23
    end do
    do ii=1,6!剪刀体6个面所在的无界面进行裁剪,外循环
        oss=0
        odd=0
        x=0.
        y=0.
        z=0.
        if(ii==1)then
            oss=6
            do i=1,oss
                odd(i)=4
                do j=1,odd(i)
                    x(i,j)=sbjx(i,j)
                    y(i,j)=sbjy(i,j)
                    z(i,j)=sbjz(i,j)
                end do
            end do
        else
            if(os==0)then
                exit
            else
                oss=os
                do i=1,oss
                    odd(i)=od(i)
                    do j=1,odd(i)
                        x(i,j)=osx(i,j)
                        y(i,j)=osy(i,j)
                        z(i,j)=osz(i,j)
                    end do
                end do
            end if
        end if
        sa=3
        sn=3
        do i=1,oss
            !判断被剪体当前所有面上所有顶点在剪刀体当前面的内侧0\面上2\外侧1
            do j=1,odd(i)
                rr=(shax(ii)-x(i,j))*(shax(ii)-jdx0)+(shay(ii)-y(i,j))*(shay(ii)-jdy0)+(shaz(ii)-z(i,j))*(shaz(ii)-jdz0)
                if(rr>0.)then!夹角小于90度,被剪顶点与剪刀重心一样在该面内
                    sa(i,j,ii)=0
                else if(rr==0.)then!夹角等于90度,被剪顶点刚好在面上
                    sa(i,j,ii)=2
                else!夹角大于90度,被剪顶点与剪刀重心在该面两侧
                    sa(i,j,ii)=1
                end if
            end do
            !判断被剪体当前所有面与剪刀体当前面是否平行
            x21=x(i,1)-x(i,2)
            y21=y(i,1)-y(i,2)
            z21=z(i,1)-z(i,2)
            x23=x(i,3)-x(i,2)
            y23=y(i,3)-y(i,2)
            z23=z(i,3)-z(i,2)
            x0=y21*z23-z21*y23  !被剪体每个面的法向量
            y0=z21*x23-x21*z23
            z0=x21*y23-y21*x23
            xx=sjdyn(ii)*z0-sjdzn(ii)*y0!剪刀体与被剪体两两面法向量叉乘
            yy=sjdzn(ii)*x0-sjdxn(ii)*z0
            zz=sjdxn(ii)*y0-sjdyn(ii)*x0
            rr=xx*xx+yy*yy+zz*zz
            if(rr==0.)then
                sn(i)=0 !两面平行
            else
                sn(i)=1 !两面不平行
            end if
        end do
        
        os=0   !裁剪后得到的面个数
        od=0   !裁剪后每个面上的点个数
        num=0  !所有面上裁剪得到的交点总数
        xj=0.
        yj=0.
        zj=0.
        do i=1,oss!被剪体的有界面,内循环
            if(sn(i)==0)then!面平行，判断其中一个点是面内上外侧即可
                if(sa(i,1,ii)==1)then  !在面外侧,全不要
                else    !在面上、内侧,全要
                    os=os+1
                    od(os)=odd(i)
                    do j=1,od(os)
                        osx(os,j)=x(i,j)
                        osy(os,j)=y(i,j)
                        osz(os,j)=z(i,j)
                    end do
                end if
            else  !两面不平行，判断面所有顶点
                k=0
                do j=1,odd(i)
                    if(sa(i,j,ii)/=0)then!所有顶点都在面外或面上，都不要
                        k=k+1
                    end if
                end do
                if(k==odd(i))then
                else
                    os=os+1
                    sum=0
                    do j=1,odd(i)
                        j1=j
                        j2=j+1
                        if(j==odd(i))then
                            j2=1
                        end if
                        if(sa(i,j1,ii)==0 .and. sa(i,j2,ii)==0)then
                            osx(os,sum+1)=x(i,j1)
                            osy(os,sum+1)=y(i,j1)
                            osz(os,sum+1)=z(i,j1)
                            osx(os,sum+2)=x(i,j2)
                            osy(os,sum+2)=y(i,j2)
                            osz(os,sum+2)=z(i,j2)
                            sum=sum+2
                        else if(sa(i,j1,ii)==0 .and. sa(i,j2,ii)==1)then
                            call interls(sjdx(ii,1),sjdy(ii,1),sjdz(ii,1),sjdx(ii,2),sjdy(ii,2),sjdz(ii,2),sjdx(ii,3),sjdy(ii,3),sjdz(ii,3),x(i,j1),y(i,j1),z(i,j1),x(i,j2),y(i,j2),z(i,j2),xx,yy,zz)
                            osx(os,sum+1)=x(i,j1)
                            osy(os,sum+1)=y(i,j1)
                            osz(os,sum+1)=z(i,j1)
                            osx(os,sum+2)=xx
                            osy(os,sum+2)=yy
                            osz(os,sum+2)=zz
                            sum=sum+2
                            xj(num+1)=xx
                            yj(num+1)=yy
                            zj(num+1)=zz
                            num=num+1
                        else if(sa(i,j1,ii)==0 .and. sa(i,j2,ii)==2)then
                            osx(os,sum+1)=x(i,j1)
                            osy(os,sum+1)=y(i,j1)
                            osz(os,sum+1)=z(i,j1)
                            osx(os,sum+2)=x(i,j2)
                            osy(os,sum+2)=y(i,j2)
                            osz(os,sum+2)=z(i,j2)
                            sum=sum+2
                            xj(num+1)=x(i,j2)
                            yj(num+1)=y(i,j2)
                            zj(num+1)=z(i,j2)
                            num=num+1
                        else if(sa(i,j1,ii)==1 .and. sa(i,j2,ii)==0)then
                            call interls(sjdx(ii,1),sjdy(ii,1),sjdz(ii,1),sjdx(ii,2),sjdy(ii,2),sjdz(ii,2),sjdx(ii,3),sjdy(ii,3),sjdz(ii,3),x(i,j1),y(i,j1),z(i,j1),x(i,j2),y(i,j2),z(i,j2),xx,yy,zz)
                            osx(os,sum+1)=xx
                            osy(os,sum+1)=yy
                            osz(os,sum+1)=zz
                            osx(os,sum+2)=x(i,j2)
                            osy(os,sum+2)=y(i,j2)
                            osz(os,sum+2)=z(i,j2)
                            sum=sum+2
                            xj(num+1)=xx
                            yj(num+1)=yy
                            zj(num+1)=zz
                            num=num+1
                        else if(sa(i,j1,ii)==1 .and. sa(i,j2,ii)==1)then
                        else if(sa(i,j1,ii)==1 .and. sa(i,j2,ii)==2)then
                        else if(sa(i,j1,ii)==2 .and. sa(i,j2,ii)==0)then
                            osx(os,sum+1)=x(i,j1)
                            osy(os,sum+1)=y(i,j1)
                            osz(os,sum+1)=z(i,j1)
                            osx(os,sum+2)=x(i,j2)
                            osy(os,sum+2)=y(i,j2)
                            osz(os,sum+2)=z(i,j2)
                            sum=sum+2
                            xj(num+1)=x(i,j1)
                            yj(num+1)=y(i,j1)
                            zj(num+1)=z(i,j1)
                            num=num+1
                        else if(sa(i,j1,ii)==2 .and. sa(i,j2,ii)==1)then
                        else if(sa(i,j1,ii)==2 .and. sa(i,j2,ii)==2)then
                            osx(os,sum+1)=x(i,j1)
                            osy(os,sum+1)=y(i,j1)
                            osz(os,sum+1)=z(i,j1)
                            osx(os,sum+2)=x(i,j2)
                            osy(os,sum+2)=y(i,j2)
                            osz(os,sum+2)=z(i,j2)
                            sum=sum+2
                            xj(num+1)=x(i,j1)
                            yj(num+1)=y(i,j1)
                            zj(num+1)=z(i,j1)
                            xj(num+2)=x(i,j2)
                            yj(num+2)=y(i,j2)
                            zj(num+2)=z(i,j2)
                            num=num+2
                        end if
                    end do
                    if(sum>1)then
                        sum0=1
                        osx0(1)=osx(os,1)
                        osy0(1)=osy(os,1)
                        osz0(1)=osz(os,1)
                        do k=2,sum
                            l=0
                            do kk=1,sum0
                                if(abs(osx0(kk)-osx(os,k))<1.d-12 .and. abs(osy0(kk)-osy(os,k))<1.d-12 .and. abs(osz0(kk)-osz(os,k))<1.d-12)then
                                    l=1
                                    exit
                                end if
                            end do
                            if(l==0)then !该点与记录里的所有点都不重合，则记录该点
                                sum0=sum0+1
                                osx0(sum0)=osx(os,k)
                                osy0(sum0)=osy(os,k)
                                osz0(sum0)=osz(os,k)
                            end if
                        end do
                        od(os)=sum0
                        osx(os,1:sum0)=osx0(1:sum0)
                        osy(os,1:sum0)=osy0(1:sum0)
                        osz(os,1:sum0)=osz0(1:sum0)
                    end if
                end if
            end if
        end do
        if(num>0)then!因裁剪出来的新界面
            i=0
            call newsec(num,xj(1:num),yj(1:num),zj(1:num),i,nnum,xjj(1:nnum),yjj(1:nnum),zjj(1:nnum))
            if(i==0)then
            else!得到新面，还要判断与某os面是否共面，一旦有就不要新面
                sum=0
                x21=xjj(1)-xjj(2)
                y21=yjj(1)-yjj(2)
                z21=zjj(1)-zjj(2)
                x23=xjj(3)-xjj(2)
                y23=yjj(3)-yjj(2)
                z23=zjj(3)-zjj(2)
                x0=y21*z23-z21*y23
                y0=z21*x23-x21*z23
                z0=x21*y23-y21*x23
                do i=1,os
                    x21=osx(i,1)-osx(i,2)
                    y21=osy(i,1)-osy(i,2)
                    z21=osz(i,1)-osz(i,2)
                    x23=osx(i,3)-osx(i,2)
                    y23=osy(i,3)-osy(i,2)
                    z23=osz(i,3)-osz(i,2)
                    xx=y21*z23-z21*y23
                    yy=z21*x23-x21*z23
                    zz=x21*y23-y21*x23
                    rrx=yy*z0-zz*y0
                    rry=zz*x0-xx*z0
                    rrz=xx*y0-yy*x0
                    rr=rrx**2+rry**2+rrz**2
                    if(rr==0.)then!两面平行
                        call shadow(xjj(1),yjj(1),zjj(1),xjj(2),yjj(2),zjj(2),xjj(3),yjj(3),zjj(3),osx(i,1),osy(i,1),osz(i,1),xx,yy,zz)
                        dd=(osx(i,1)-xx)**2+(osy(i,1)-yy)**2+(osz(i,1)-zz)**2
                        if(dd==0.)then!两面重合
                            sum=1
                            exit
                        end if
                    end if
                end do
                if(sum==0)then
                    os=os+1
                    od(os)=nnum
                    do j=1,od(os)
                        osx(os,j)=xjj(j)
                        osy(os,j)=yjj(j)
                        osz(os,j)=zjj(j)
                    end do
                end if
            end if
        end if
    end do
    !求多面体体积
    v=0.
    if(os>=4)then
        !求多面体的体中心，要找出该体的所有不重复顶点
        do j=1,od(1)
            zdx(j)=osx(1,j)
            zdy(j)=osy(1,j)
            zdz(j)=osz(1,j)
        end do
        nnum=od(1)
        do i=2,os
            do j=1,od(i)
                l=0
                do k=1,nnum
                    if(abs(osx(i,j)-zdx(k))<1.d-12 .and. abs(osy(i,j)-zdy(k))<1.d-12 .and. abs(osz(i,j)-zdz(k))<1.d-12)then
                        l=1
                        exit
                    end if
                end do
                if(l==0)then !该点与记录里的所有点都不重合，则记录该点
                    nnum=nnum+1
                    zdx(nnum)=osx(i,j)
                    zdy(nnum)=osy(i,j)
                    zdz(nnum)=osz(i,j)
                end if
            end do
        end do
        x0=0.
        y0=0.
        z0=0.
        do i=1,nnum
            x0=x0+zdx(i)
            y0=y0+zdy(i)
            z0=z0+zdz(i)
        end do
        x0=x0/dble(nnum)
        y0=y0/dble(nnum)
        z0=z0/dble(nnum)
        do i=1,os
            call shadow(osx(i,1),osy(i,1),osz(i,1),osx(i,2),osy(i,2),osz(i,2),osx(i,3),osy(i,3),osz(i,3),x0,y0,z0,xx,yy,zz)
            hh=sqrt((xx-x0)**2+(yy-y0)**2+(zz-z0)**2)
            dd=0.
            do j=3,od(i)
                x21=osx(i,j-1)-osx(i,1)
                y21=osy(i,j-1)-osy(i,1)
                z21=osz(i,j-1)-osz(i,1)
                x23=osx(i,j)-osx(i,1)
                y23=osy(i,j)-osy(i,1)
                z23=osz(i,j)-osz(i,1)
                xx=y21*z23-z21*y23
                yy=z21*x23-x21*z23
                zz=x21*y23-y21*x23
                rr=sqrt(xx*xx+yy*yy+zz*zz)
                dd=dd+rr
            end do
            v=v+hh*dd/6.
        end do
    end if
    end subroutine clipv
    
    subroutine newsec(num,xj,yj,zj,ii,nnum,x,y,z)
    implicit none
    integer :: num
    real(8) :: xj(num),yj(num),zj(num)
    integer :: ii,nnum
    real(8) :: x(nnum),y(nnum),z(nnum)
    real(8) :: xjj(10),yjj(10),zjj(10),xj0(10),yj0(10),zj0(10),xj1(10),yj1(10),zj1(10),xj2(10),yj2(10),zj2(10),tha(10),xn(10),yn(10),zn(10)
    real(8) :: xj01(10),yj01(10),zj01(10),xj02(10),yj02(10),zj02(10),xj03(10),yj03(10),zj03(10),rr0(10),rr01(10),rr02(10),rr03(10)
    integer :: i,j,k,l,jj,kk,k1,k2,k3
    real(8) :: x0,y0,z0,x21,y21,z21,rr21,x23,y23,z23,rr23,temp,rr
    
    nnum=1
    xjj(1)=xj(1)
    yjj(1)=yj(1)
    zjj(1)=zj(1)
    do j=2,num
        l=0
        do k=1,nnum
            if(abs(xj(j)-xjj(k))<1.d-12 .and. abs(yj(j)-yjj(k))<1.d-12 .and. abs(zj(j)-zjj(k))<1.d-12)then
                l=1
                exit
            end if
        end do
        if(l==0)then !该点与记录里的所有点都不重合，则记录该点
            nnum=nnum+1
            xjj(nnum)=xj(j)
            yjj(nnum)=yj(j)
            zjj(nnum)=zj(j)
        end if
    end do
    
    if(nnum<3)then
        ii=0
        x=0.
        y=0.
        z=0.
    else
        ii=1
        !构成新界面，对点排序
        x0=0.
        y0=0.
        z0=0.
        do j=1,nnum
            x0=x0+xjj(j)
            y0=y0+yjj(j)
            z0=z0+zjj(j)
        end do
        x0=x0/dble(nnum)!多边形中心点
        y0=y0/dble(nnum)
        z0=z0/dble(nnum)
        !先将所有点关于是否在OA一侧进行分类（A为其中一点）
        kk=0 !在OA直线上
        do j=2,nnum
            xn(j)=(y0-yjj(j))*(zjj(1)-zjj(j))-(z0-zjj(j))*(yjj(1)-yjj(j))
            yn(j)=(z0-zjj(j))*(xjj(1)-xjj(j))-(x0-xjj(j))*(zjj(1)-zjj(j))
            zn(j)=(x0-xjj(j))*(yjj(1)-yjj(j))-(y0-yjj(j))*(xjj(1)-xjj(j))
            rr=xn(j)**2+yn(j)**2+zn(j)**2
            if(rr==0.)then
                kk=kk+1
                xj0(kk)=xjj(j)
                yj0(kk)=yjj(j)
                zj0(kk)=zjj(j)
            end if
        end do
        k=1  !假设2点在k侧，判断所有点与2点是否在一侧
        l=0  !若与2点不在一侧，则在l侧
        xj1(1)=xjj(2)
        yj1(1)=yjj(2)
        zj1(1)=zjj(2)
        do j=3,nnum
            rr=xn(2)*xn(j)+yn(2)*yn(j)+zn(2)*zn(j)
            if(rr>0.)then!该点与2点在OA一侧
                k=k+1
                xj1(k)=xjj(j)
                yj1(k)=yjj(j)
                zj1(k)=zjj(j)
            else if(rr<0.)then!该点与2点在OA另一侧
                l=l+1
                xj2(l)=xjj(j)
                yj2(l)=yjj(j)
                zj2(l)=zjj(j)
            end if
        end do
        !每一侧按夹角进行排序（某一侧按0到pai排序，另一侧接着按pai到0排序）
        x21=xjj(1)-x0 !基准边
        y21=yjj(1)-y0
        z21=zjj(1)-z0
        rr21=sqrt(x21*x21+y21*y21+z21*z21)
        do j=1,k
            x23=xj1(j)-x0
            y23=yj1(j)-y0
            z23=zj1(j)-z0
            rr23=sqrt(x23*x23+y23*y23+z23*z23)
            rr=x21*x23+y21*y23+z21*z23  !向量点积
            tha(j)=acos(rr/rr21/rr23) !两向量夹角
        end do
        do j=1,k-1
            do jj=j+1,k
                if(tha(j)>tha(jj))then
                    temp=tha(j)
                    tha(j)=tha(jj)
                    tha(jj)=temp
                    temp=xj1(j)
                    xj1(j)=xj1(jj)
                    xj1(jj)=temp
                    temp=yj1(j)
                    yj1(j)=yj1(jj)
                    yj1(jj)=temp
                    temp=zj1(j)
                    zj1(j)=zj1(jj)
                    zj1(jj)=temp
                end if
            end do
        end do
        do j=1,l
            x23=xj2(j)-x0
            y23=yj2(j)-y0
            z23=zj2(j)-z0
            rr23=sqrt(x23*x23+y23*y23+z23*z23)
            rr=x21*x23+y21*y23+z21*z23  !向量点积
            tha(j)=acos(rr/rr21/rr23) !两向量夹角
        end do
        do j=1,l-1
            do jj=j+1,l
                if(tha(j)<tha(jj))then
                    temp=tha(j)
                    tha(j)=tha(jj)
                    tha(jj)=temp
                    temp=xj2(j)
                    xj2(j)=xj2(jj)
                    xj2(jj)=temp
                    temp=yj2(j)
                    yj2(j)=yj2(jj)
                    yj2(jj)=temp
                    temp=zj2(j)
                    zj2(j)=zj2(jj)
                    zj2(jj)=temp
                end if
            end do
        end do
        !OA线上点排序
        k1=0
        k2=0
        k3=0
        do j=1,kk
            x23=xj0(j)-x0
            y23=yj0(j)-y0
            z23=zj0(j)-z0
            rr0(j)=sqrt(x23*x23+y23*y23+z23*z23)
            rr=x23*x21+y23*y21+z23*z21
            if(rr<0.)then!在O左侧
                k1=k1+1
                xj01(k1)=xj0(j)
                yj01(k1)=yj0(j)
                zj01(k1)=zj0(j)
                rr01(k1)=rr0(j)
            else
                if(rr0(j)<rr21)then!在OA内
                    k2=k2+1
                    xj02(k2)=xj0(j)
                    yj02(k2)=yj0(j)
                    zj02(k2)=zj0(j)
                    rr02(k2)=rr0(j)
                else               !在A右侧
                    k3=k3+1
                    xj03(k3)=xj0(j)
                    yj03(k2)=yj0(j)
                    zj03(k3)=zj0(j)
                    rr03(k3)=rr0(j)
                end if
            end if
        end do
        do j=1,k1-1
            do jj=j+1,k1
                if(rr01(j)>rr01(jj))then
                    temp=rr01(j)
                    rr01(j)=rr01(jj)
                    rr01(jj)=temp
                    temp=xj01(j)
                    xj01(j)=xj01(jj)
                    xj01(jj)=temp
                    temp=yj01(j)
                    yj01(j)=yj01(jj)
                    yj01(jj)=temp
                    temp=zj01(j)
                    zj01(j)=zj01(jj)
                    zj01(jj)=temp
                end if
            end do
        end do
        do j=1,k2-1
            do jj=j+1,k2
                if(rr02(j)>rr02(jj))then
                    temp=rr02(j)
                    rr02(j)=rr02(jj)
                    rr02(jj)=temp
                    temp=xj02(j)
                    xj02(j)=xj02(jj)
                    xj02(jj)=temp
                    temp=yj02(j)
                    yj02(j)=yj02(jj)
                    yj02(jj)=temp
                    temp=zj02(j)
                    zj02(j)=zj02(jj)
                    zj02(jj)=temp
                end if
            end do
        end do
        do j=1,k3-1
            do jj=j+1,k3
                if(rr03(j)<rr03(jj))then
                    temp=rr03(j)
                    rr03(j)=rr03(jj)
                    rr03(jj)=temp
                    temp=xj03(j)
                    xj03(j)=xj03(jj)
                    xj03(jj)=temp
                    temp=yj03(j)
                    yj03(j)=yj03(jj)
                    yj03(jj)=temp
                    temp=zj03(j)
                    zj03(j)=zj03(jj)
                    zj03(jj)=temp
                end if
            end do
        end do
        !重组
        x(1)=xjj(1)
        y(1)=yjj(1)
        z(1)=zjj(1)
        do j=1,k
            x(j+1)=xj1(j)
            y(j+1)=yj1(j)
            z(j+1)=zj1(j)
        end do
        do j=1,k1
            x(j+k+1)=xj01(j)
            y(j+k+1)=yj01(j)
            z(j+k+1)=zj01(j)
        end do
        do j=1,l
            x(j+k+k1+1)=xj2(j)
            y(j+k+k1+1)=yj2(j)
            z(j+k+k1+1)=zj2(j)
        end do
        do j=1,k2
            x(j+k+k1+l+1)=xj02(j)
            y(j+k+k1+l+1)=yj02(j)
            z(j+k+k1+l+1)=zj02(j)
        end do
        do j=1,k3
            x(j+k+k1+l+k2+1)=xj03(j)
            y(j+k+k1+l+k2+1)=yj03(j)
            z(j+k+k1+l+k2+1)=zj03(j)
        end do
    end if
    end subroutine newsec
    
    subroutine shadow(x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0,x00,y00,z00)  !点到面的垂足
    implicit none
    real(8),intent(in) :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0
    real(8),intent(out) :: x00,y00,z00
    real(8) :: x21,y21,z21,x23,y23,z23,xx,yy,zz,rr,rc
    
    x21=x1-x2
    y21=y1-y2
    z21=z1-z2
    x23=x3-x2
    y23=y3-y2
    z23=z3-z2
    xx=y21*z23-z21*y23  !法向量
    yy=z21*x23-x21*z23
    zz=x21*y23-y21*x23
    rr=xx*xx+yy*yy+zz*zz
    if(rr==0.)then
        rr=1.d-20
    end if
    rc=((x1-x0)*xx+(y1-y0)*yy+(z1-z0)*zz)/rr
    x00=x0+xx*rc!同一方向上知道长度,确定坐标
    y00=y0+yy*rc
    z00=z0+zz*rc
    end subroutine shadow
    
    subroutine interls(sx1,sy1,sz1,sx2,sy2,sz2,sx3,sy3,sz3,x1,y1,z1,x2,y2,z2,x0,y0,z0)  !线段与面的交点
    implicit none
    real(8),intent(in) :: sx1,sy1,sz1,sx2,sy2,sz2,sx3,sy3,sz3,x1,y1,z1,x2,y2,z2
    real(8),intent(out) :: x0,y0,z0
    real(8) :: x10,y10,z10,x20,y20,z20,rr1,rr2,xx,yy,zz,rc

    call shadow(sx1,sy1,sz1,sx2,sy2,sz2,sx3,sy3,sz3,x1,y1,z1,x10,y10,z10)!垂足
    call shadow(sx1,sy1,sz1,sx2,sy2,sz2,sx3,sy3,sz3,x2,y2,z2,x20,y20,z20)
    rr1=sqrt((x10-x1)**2+(y10-y1)**2+(z10-z1)**2)
    rr2=sqrt((x20-x2)**2+(y20-y2)**2+(z20-z2)**2)
    rc=rr1/(rr1+rr2)
    xx=x20-x10
    yy=y20-y10
    zz=z20-z10
    x0=x10+xx*rc  !求出交点
    y0=y10+yy*rc
    z0=z10+zz*rc
    end subroutine interls
    
subroutine tsd!二阶向后差分中的前两时层值
    use global
    implicit none
    
    if(m==1 .and. n==2)then
        aa3=1.
        aa2=-1.
        aa1=0
    else
        aa3=1.5
        aa2=-2.
        aa1=0.5
    end if
    do k=1,nz
        do j=1,ny
            do i=1,nx
                ts1(i,j,k)=aa2*q21(i,j,k,2)+aa1*q21(i,j,k,1)
                ts2(i,j,k)=aa2*q22(i,j,k,2)+aa1*q22(i,j,k,1)
                ts3(i,j,k)=aa2*q23(i,j,k,2)+aa1*q23(i,j,k,1)
                ts4(i,j,k)=aa2*q24(i,j,k,2)+aa1*q24(i,j,k,1)
                ts5(i,j,k)=aa2*q25(i,j,k,2)+aa1*q25(i,j,k,1)
                ts6(i,j,k)=aa2*q26(i,j,k,2)+aa1*q26(i,j,k,1)
            end do
        end do
    end do
end subroutine tsd

subroutine march !3步R-K推进
    use global
    implicit none
    
    q01(1:nx,1:ny,1:nz)=q11(1:nx,1:ny,1:nz)
    q02(1:nx,1:ny,1:nz)=q12(1:nx,1:ny,1:nz)
    q03(1:nx,1:ny,1:nz)=q13(1:nx,1:ny,1:nz)
    q04(1:nx,1:ny,1:nz)=q14(1:nx,1:ny,1:nz)
    q05(1:nx,1:ny,1:nz)=q15(1:nx,1:ny,1:nz)
    q06(1:nx,1:ny,1:nz)=q16(1:nx,1:ny,1:nz)

    timl=0.6d0
    call ppp
    call bc
    call step
    call ddd
    call qqq
    call qqqv
    call pred(1)
    
    timl=0.6d0
    call ppp
    call bc
    call qqq
    call pred(0)
    
    timl=1.d0
    call ppp
    call bc
    call qqq
    call pred(1)
    
    q01(1:nx,1:ny,1:nz)=q11(1:nx,1:ny,1:nz)
    q02(1:nx,1:ny,1:nz)=q12(1:nx,1:ny,1:nz)
    q03(1:nx,1:ny,1:nz)=q13(1:nx,1:ny,1:nz)
    q04(1:nx,1:ny,1:nz)=q14(1:nx,1:ny,1:nz)
    q05(1:nx,1:ny,1:nz)=q15(1:nx,1:ny,1:nz)
    q06(1:nx,1:ny,1:nz)=q16(1:nx,1:ny,1:nz)

    timl=0.125d0
    call ppp
    call bc
    call step
    call ddd
    call qqq
    call qqqv
    call pred(0)
    end subroutine march

subroutine ppp  !由守恒变量求原始变量
    use global
    implicit none
    
    do k=1,nz
        do j=1,ny
            do i=1,nx
                dim=Q11(i,j,k)
                pvx(i,j,k)=Q12(i,j,k)/dim
                pvy(i,j,k)=Q13(i,j,k)/dim
                pvz(i,j,k)=Q14(i,j,k)/dim
                vx=pvx(i,j,k)
                vy=pvy(i,j,k)
                vz=pvz(i,j,k)
                
                y1=yy0(i,j,k)
                z1=zz0(i,j,k)
                rr=sqrt(y1*y1+z1*z1)
                sir=z1/rr
                cor=y1/rr
                vth(i,j,k)=vz*cor-vy*sir    !叶高y速度
                vre(i,j,k)=vz*sir+vy*cor    !周向z速度
                
                qq2=vx*vx+vy*vy+vz*vz
                p(i,j,k)=0.4d0*(Q15(i,j,k)-0.5d0*dim*qq2)
                t(i,j,k)=p(i,j,k)/(dim*rg)
                cvl=cvl0*(t(i,j,k)/t0)**1.5*(t0+ts)/(t(i,j,k)+ts)
                q16(i,j,k)=max(q16(i,j,k),1.D-4*cvl)
            end do
        end do
    end do
    end subroutine ppp
    
subroutine bc
    use global
    implicit none
    real(8) ::sxn,syn,szn,ds,deltp,rrhoc
    real(8) ::uabs,unorm,rinv,c02,dis,cb,cosa,hb,cc02
    real(8) ::sss,sr,sv,sd,s1,ve,v2,dr,x1,y2,z2,rr22,rp1,rp2
    real(8) ::ca,qz1,qz2,qq,t1,cvlt,cvu,tem,tur1,tur2,tur3,fv1,temp
    real(8),allocatable ::hr(:),hv(:),hd(:),sq11(:,:,:,:),spvx(:,:,:,:),spvy(:,:,:,:),spvz(:,:,:,:),sp(:,:,:,:),sq16(:,:,:,:)
    integer :: num,iii,cycle_by
    !****开始通信，待通信结束后再把通信后的量转换为守恒量
    !****x方向动静叶内交界面传递
    if(slidm/=0)then
        allocate(sq11(2,1:ny,1:nz,3))
        allocate(spvx(2,1:ny,1:nz,3))
        allocate(spvy(2,1:ny,1:nz,3))
        allocate(spvz(2,1:ny,1:nz,3))
        allocate(sp(2,1:ny,1:nz,3))
        allocate(sq16(2,1:ny,1:nz,3))
        if(rm==1)then
            i=nx-1
        else if(rm==2)then
            i=1
        end if
        do l=1,ssum
            num=cli(l)
            call MPI_SENDRECV(q11(i:i+1,1:ny,bkk(l):ekk(l)),2*ny*(ekk(l)-bkk(l)+1),MPI_DOUBLE_PRECISION,num,115,sq11(1:2,1:ny,bk(l):ek(l),l),2*ny*(ek(l)-bk(l)+1),MPI_DOUBLE_PRECISION,num,115,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(pvx(i:i+1,1:ny,bkk(l):ekk(l)),2*ny*(ekk(l)-bkk(l)+1),MPI_DOUBLE_PRECISION,num,116,spvx(1:2,1:ny,bk(l):ek(l),l),2*ny*(ek(l)-bk(l)+1),MPI_DOUBLE_PRECISION,num,116,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(pvy(i:i+1,1:ny,bkk(l):ekk(l)),2*ny*(ekk(l)-bkk(l)+1),MPI_DOUBLE_PRECISION,num,117,spvy(1:2,1:ny,bk(l):ek(l),l),2*ny*(ek(l)-bk(l)+1),MPI_DOUBLE_PRECISION,num,117,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(pvz(i:i+1,1:ny,bkk(l):ekk(l)),2*ny*(ekk(l)-bkk(l)+1),MPI_DOUBLE_PRECISION,num,118,spvz(1:2,1:ny,bk(l):ek(l),l),2*ny*(ek(l)-bk(l)+1),MPI_DOUBLE_PRECISION,num,118,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(p(i:i+1,1:ny,bkk(l):ekk(l)),2*ny*(ekk(l)-bkk(l)+1),MPI_DOUBLE_PRECISION,num,119,sp(1:2,1:ny,bk(l):ek(l),l),2*ny*(ek(l)-bk(l)+1),MPI_DOUBLE_PRECISION,num,119,MPI_COMM_WORLD,status,ierr)
            call MPI_SENDRECV(q16(i:i+1,1:ny,bkk(l):ekk(l)),2*ny*(ekk(l)-bkk(l)+1),MPI_DOUBLE_PRECISION,num,120,sq16(1:2,1:ny,bk(l):ek(l),l),2*ny*(ek(l)-bk(l)+1),MPI_DOUBLE_PRECISION,num,120,MPI_COMM_WORLD,status,ierr)
        end do
    end if
    !***x方向叶片内进程传递信息
    if(xll>0)then
        if(xl==1)then
            j=myid-1
        else
            if(yl==0)then
                j=myid-xln(xlln-1)-1
            else
                j=myid-xln(xlln)-1
            end if
        end if
        call MPI_SENDRECV(q11(1,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,570,q11(0,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,560,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvx(1,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,571,pvx(0,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,561,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvy(1,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,572,pvy(0,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,562,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvz(1,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,573,pvz(0,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,563,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(p(1,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,574,p(0,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,564,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(q16(1,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,575,q16(0,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,565,MPI_COMM_WORLD,status,ierr)
    end if
    if(xll<xln(0)+xln(1)+xln(2)-1)then
        if(xl==0)then
            if(slidm==2)then
                if(yl==0)then
                    j=myid+xln(xlln)+1
                else
                    j=myid+xln(xlln+1)+1
                end if
            else
                j=myid+1
            end if
        else
            if(yl==0)then
                j=myid+xln(xlln)+1
            else
                j=myid+xln(xlln+1)+1
            end if
        end if
        call MPI_SENDRECV(q11(nx,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,560,q11(nx+1,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,570,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvx(nx,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,561,pvx(nx+1,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,571,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvy(nx,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,562,pvy(nx+1,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,572,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvz(nx,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,563,pvz(nx+1,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,573,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(p(nx,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,564,p(nx+1,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,574,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(q16(nx,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,565,q16(nx+1,1:ny,1:nz),ny*nz,MPI_DOUBLE_PRECISION,j,575,MPI_COMM_WORLD,status,ierr)
    end if
    !***y方向叶片内传递信息
    if(yl<yln-1)then
        call MPI_SENDRECV(q11(1:nx,ny,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),20,q11(1:nx,ny+1,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),30,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvx(1:nx,ny,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),21,pvx(1:nx,ny+1,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),31,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvy(1:nx,ny,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),22,pvy(1:nx,ny+1,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),32,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvz(1:nx,ny,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),23,pvz(1:nx,ny+1,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),33,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(p(1:nx,ny,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),24,p(1:nx,ny+1,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),34,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(q16(1:nx,ny,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),25,q16(1:nx,ny+1,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid+xln(xlln),35,MPI_COMM_WORLD,status,ierr)
    end if
    if(yl>0)then
        call MPI_SENDRECV(q11(1:nx,1,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),30,q11(1:nx,0,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),20,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvx(1:nx,1,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),31,pvx(1:nx,0,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),21,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvy(1:nx,1,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),32,pvy(1:nx,0,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),22,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvz(1:nx,1,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),33,pvz(1:nx,0,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),23,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(p(1:nx,1,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),34,p(1:nx,0,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),24,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(q16(1:nx,1,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),35,q16(1:nx,0,1:nz),nx*nz,MPI_DOUBLE_PRECISION,myid-xln(xlln),25,MPI_COMM_WORLD,status,ierr)
    end if
    !***********z方向各通道间
    if(lm==0)then
        myidl=myid+(lb(rm)-1)*numpp
        myidr=myid+numpp
    else if(lm==lb(rm)-1)then
        myidl=myid-numpp
        myidr=myid-(lb(rm)-1)*numpp
    else
        myidl=myid-numpp
        myidr=myid+numpp
    end if
    if(xlln==1)then
        cycle_by=jt+1
    else
        cycle_by=1
    end if
    if(cycle_by<=ny)then
        call MPI_SENDRECV(q11(1:nx,cycle_by:ny,1),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidr,120,q11(1:nx,cycle_by:ny,nz+1),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidl,120,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvx(1:nx,cycle_by:ny,1),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidr,121,pvx(1:nx,cycle_by:ny,nz+1),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidl,121,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvy(1:nx,cycle_by:ny,1),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidr,122,pvy(1:nx,cycle_by:ny,nz+1),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidl,122,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvz(1:nx,cycle_by:ny,1),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidr,123,pvz(1:nx,cycle_by:ny,nz+1),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidl,123,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(p(1:nx,cycle_by:ny,1),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidr,124,p(1:nx,cycle_by:ny,nz+1),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidl,124,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(q16(1:nx,cycle_by:ny,1),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidr,125,q16(1:nx,cycle_by:ny,nz+1),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidl,125,MPI_COMM_WORLD,status,ierr)
        
        call MPI_SENDRECV(q11(1:nx,cycle_by:ny,nz),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidl,130,q11(1:nx,cycle_by:ny,0),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidr,130,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvx(1:nx,cycle_by:ny,nz),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidl,131,pvx(1:nx,cycle_by:ny,0),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidr,131,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvy(1:nx,cycle_by:ny,nz),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidl,132,pvy(1:nx,cycle_by:ny,0),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidr,132,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(pvz(1:nx,cycle_by:ny,nz),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidl,133,pvz(1:nx,cycle_by:ny,0),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidr,133,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(p(1:nx,cycle_by:ny,nz),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidl,134,p(1:nx,cycle_by:ny,0),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidr,134,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(q16(1:nx,cycle_by:ny,nz),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidl,135,q16(1:nx,cycle_by:ny,0),nx*(ny-cycle_by+1),MPI_DOUBLE_PRECISION,myidr,135,MPI_COMM_WORLD,status,ierr)
    end if
    !***************初级进口涡粘度turi分布、入口第一进程的进口边界***************
    if(rm==1 .and. xll==0)then
        do k=1,nz
            do j=1,ny
                dim=Q11(1,j,k)
                vx=Q12(1,j,k)/dim
                vy=Q13(1,j,k)/dim
                vz=Q14(1,j,k)/dim
                qq2=vx*vx+vy*vy+vz*vz
                en=Q15(1,j,k)
                pp=0.4d0*(en-0.5d0*dim*qq2)
                t1=pp/(dim*rg)
                cvlt=cvl0*(t1/t0)**1.5*(t0+ts)/(t1+ts)
                cvl=cvlt/dim
                cvu=c2*cvl          !重点是c2的给值
                tur1=1.D-4
                tur2=1d-6
                tur3=1d-6
                do while(abs((tur1-tur3)/tur3)>1d-6)
                    tem=tur1/cvl
                    fv1=1.d0/(1.d0+(cv1/tem)**3)
                    tur2=cvu/fv1
                    tur3=tur1
                    tur1=tur2
                end do
                
                ds=sqrt(s2x(1,j,k)*s2x(1,j,k)+s2y(1,j,k)*s2y(1,j,k)+s2z(1,j,k)*s2z(1,j,k))
                sxn=s2x(1,j,k)/ds
                syn=s2y(1,j,k)/ds
                szn=s2z(1,j,k)/ds
                vx=pvx(1,j,k)
                vy=pvy(1,j,k)
                vz=pvz(1,j,k)
                uabs=sqrt(vx*vx+vy*vy+vz*vz)
                unorm = vx*sxn + vy*syn +vz*szn
                a=sqrt(1.4d0*p(1,j,k)/q11(1,j,k))
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
                hb    = ht*cc02
            
                t(0,j,k)=hb/cp
                p(0,j,k)=pt*(cc02**3.5)
                q11(0,j,k)=p(0,j,k)/(t(0,j,k)*rg)
                uabs =2.d0*(ht-hb)
                q15(0,j,k) =2.5d0*p(0,j,k)+0.5d0*q11(0,j,k)*uabs
                pvx(0,j,k)=sqrt(uabs)*betax(j,k)               !入口速度垂直入口
                pvy(0,j,k)=sqrt(uabs)*betay(j,k)
                pvz(0,j,k)=sqrt(uabs)*betaz(j,k)
                q12(0,j,k)=q11(0,j,k)*pvx(0,j,k)
                q13(0,j,k)=q11(0,j,k)*pvy(0,j,k)
                q14(0,j,k)=q11(0,j,k)*pvz(0,j,k)
                q16(0,j,k)=tur2*q11(0,j,k)
            end do
        end do
    end if
    !***************末级出口静压peb分布、出口进程的出口边界***************
    if(rm==2 .and. xll==xln(0)+xln(1)+xln(2)-1)then
          allocate(hr(0:ny))
          allocate(hv(0:ny))
           allocate(hd(0:ny))
        do j=1,ny
            sss=0.d0
            sr=0.d0
            sv=0.d0
            sd=0.d0
            do k=1,nz
                s1=sqrt(s2x(nx+1,j,k)*s2x(nx+1,j,k)+s2y(nx+1,j,k)*s2y(nx+1,j,k)+s2z(nx+1,j,k)*s2z(nx+1,j,k))
                y1=yy02(nx+1,j,k)
                z1=zz02(nx+1,j,k)
                rr=sqrt(y1*y1+z1*z1)
                dim=q11(nx,j,k)
                ve=vth(nx,j,k)
                sss=sss+s1
                sr=sr+rr*s1
                sd=sd+dim*s1
                sv=sv+ve*s1
            end do
            hr(j)=sr/sss
            hv(j)=sv/sss
            hd(j)=sd/sss
        end do
        if(yl==0)then
               peb(1)=pb1
            do j=2,ny
                dim=0.5d0*(hd(j-1)+hd(j))
                v2=0.5d0*(hv(j-1)+hv(j))
                dr=hr(j)-hr(j-1)
                rr=0.5d0*(hr(j)+hr(j-1))
                peb(j)=peb(j-1)+dim*v2*v2/rr*dr
            end do
            call MPI_SEND(peb(ny),1,MPI_DOUBLE_PRECISION,myid+xln(2),3,MPI_COMM_WORLD,ierr)
            call MPI_SEND(hd(ny),1,MPI_DOUBLE_PRECISION,myid+xln(2),4,MPI_COMM_WORLD,ierr)
            call MPI_SEND(hv(ny),1,MPI_DOUBLE_PRECISION,myid+xln(2),5,MPI_COMM_WORLD,ierr)
            call MPI_SEND(hr(ny),1,MPI_DOUBLE_PRECISION,myid+xln(2),6,MPI_COMM_WORLD,ierr)
        else if(yl==1)then
            call MPI_RECV(peb(0),1,MPI_DOUBLE_PRECISION,myid-xln(2),3,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(hd(0),1,MPI_DOUBLE_PRECISION,myid-xln(2),4,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(hv(0),1,MPI_DOUBLE_PRECISION,myid-xln(2),5,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(hr(0),1,MPI_DOUBLE_PRECISION,myid-xln(2),6,MPI_COMM_WORLD,status,ierr)
            do j=1,ny
                dim=0.5d0*(hd(j-1)+hd(j))
                v2=0.5d0*(hv(j-1)+hv(j))
                dr=hr(j)-hr(j-1)
                rr=0.5d0*(hr(j)+hr(j-1))
                peb(j)=peb(j-1)+dim*v2*v2/rr*dr
            end do
        end if
        deallocate(hr)
        deallocate(hv)
        deallocate(hd)
        
        do k=1,nz
            do j=1,ny
                ds=sqrt(s2x(nx+1,j,k)*s2x(nx+1,j,k)+s2y(nx+1,j,k)*s2y(nx+1,j,k)+s2z(nx+1,j,k)*s2z(nx+1,j,k))
                sxn=-s2x(nx+1,j,k)/ds
                syn=-s2y(nx+1,j,k)/ds
                szn=-s2z(nx+1,j,k)/ds
                a=sqrt(1.4d0*p(nx,j,k)/q11(nx,j,k))
                rrhoc=1.d0/(q11(nx,j,k)*a)
                deltp=p(nx,j,k)-peb(j)
            
                p(nx+1,j,k)=peb(j)
                q11(nx+1,j,k) =q11(nx,j,k)-deltp/(a*a)
                pvx(nx+1,j,k) =pvx(nx,j,k)+sxn*deltp*rrhoc
                pvy(nx+1,j,k) =pvy(nx,j,k)+syn*deltp*rrhoc
                pvz(nx+1,j,k) =pvz(nx,j,k)+szn*deltp*rrhoc
                q12(nx+1,j,k)=q11(nx+1,j,k)*pvx(nx+1,j,k)
                q13(nx+1,j,k)=q11(nx+1,j,k)*pvy(nx+1,j,k)
                q14(nx+1,j,k)=q11(nx+1,j,k)*pvz(nx+1,j,k)
                qq2=pvx(nx+1,j,k)*pvx(nx+1,j,k)+pvy(nx+1,j,k)*pvy(nx+1,j,k)+pvz(nx+1,j,k)*pvz(nx+1,j,k)
                q15(nx+1,j,k) =2.5d0*p(nx+1,j,k)+0.5d0*q11(nx+1,j,k)*qq2
                q16(nx+1,j,k) =q16(nx,j,k)
                t(nx+1,j,k)=p(nx+1,j,k)/(q11(nx+1,j,k)*rg)
            end do
        end do
    end if
    !***************上下边界y***************
    if(yl==0)then
        do i=1,nx
            do k=1,nz
                !*****叶根位置****动叶轮鼓一直在转
                p(i,0,k)  = p(i,1,k)
                t(i,0,k)  = t(i,1,k)
                q11(i,0,k)= q11(i,1,k)
                pvx(i,0,k)=-pvx(i,1,k)
                pvy(i,0,k)=-2.d0*rpm(rm)*zz03(i,1,k)-pvy(i,1,k)
                pvz(i,0,k)= 2.d0*rpm(rm)*yy03(i,1,k)-pvz(i,1,k)
                q12(i,0,k)=q11(i,0,k)*pvx(i,0,k)
                q13(i,0,k)=q11(i,0,k)*pvy(i,0,k)
                q14(i,0,k)=q11(i,0,k)*pvz(i,0,k)
                qq2=pvx(i,0,k)*pvx(i,0,k)+pvy(i,0,k)*pvy(i,0,k)+pvz(i,0,k)*pvz(i,0,k)
                q15(i,0,k)=2.5d0*p(i,0,k)+0.5d0*q11(i,0,k)*qq2
                q16(i,0,k)=-q16(i,1,k)
            end do
        end do
    end if
    if(yl==yln-1)then
        do i=1,nx
            do k=1,nz
                !*****叶顶位置****
                p(i,ny+1,k)  = p(i,ny,k)
                t(i,ny+1,k)  = t(i,ny,k)
                q11(i,ny+1,k)=q11(i,ny,k)
                pvx(i,ny+1,k)=-pvx(i,ny,k)
                pvy(i,ny+1,k)=-pvy(i,ny,k)
                pvz(i,ny+1,k)=-pvz(i,ny,k)
                q12(i,ny+1,k)=q11(i,ny+1,k)*pvx(i,ny+1,k)
                q13(i,ny+1,k)=q11(i,ny+1,k)*pvy(i,ny+1,k)
                q14(i,ny+1,k)=q11(i,ny+1,k)*pvz(i,ny+1,k)
                qq2=pvx(i,ny+1,k)*pvx(i,ny+1,k)+pvy(i,ny+1,k)*pvy(i,ny+1,k)+pvz(i,ny+1,k)*pvz(i,ny+1,k)
                q15(i,ny+1,k)=2.5d0*p(i,ny+1,k)+0.5d0*q11(i,ny+1,k)*qq2
                q16(i,ny+1,k)=-q16(i,ny,k)
            end do
        end do
    end if
    !***************前后固壁边界z的计算***************
    if(xlln==1)then
        do j=jb,jt
            do i=1,nx
                p(i,j,0)=p(i,j,1)
                t(i,j,0)=t(i,j,1)
                q11(i,j,0)=q11(i,j,1)
                pvx(i,j,0)=-pvx(i,j,1)
                pvy(i,j,0)=-2.d0*rpm(rm)*zz01(i,j,1)-pvy(i,j,1)
                pvz(i,j,0)= 2.d0*rpm(rm)*yy01(i,j,1)-pvz(i,j,1)
                q12(i,j,0)=q11(i,j,0)*pvx(i,j,0)
                q13(i,j,0)=q11(i,j,0)*pvy(i,j,0)
                q14(i,j,0)=q11(i,j,0)*pvz(i,j,0)
                qq2=pvx(i,j,0)*pvx(i,j,0)+pvy(i,j,0)*pvy(i,j,0)+pvz(i,j,0)*pvz(i,j,0)
                q15(i,j,0)=2.5d0*p(i,j,0)+0.5d0*q11(i,j,0)*qq2
                q16(i,j,0)=-q16(i,j,1)
                
                p(i,j,nz+1)=p(i,j,nz)
                t(i,j,nz+1)=t(i,j,nz)
                q11(i,j,nz+1)=q11(i,j,nz)
                pvx(i,j,nz+1)=-pvx(i,j,nz)
                pvy(i,j,nz+1)=-2.d0*rpm(rm)*zz01(i,j,nz+1)-pvy(i,j,nz)
                pvz(i,j,nz+1)= 2.d0*rpm(rm)*yy01(i,j,nz+1)-pvz(i,j,nz)
                q12(i,j,nz+1)=q11(i,j,nz+1)*pvx(i,j,nz+1)
                q13(i,j,nz+1)=q11(i,j,nz+1)*pvy(i,j,nz+1)
                q14(i,j,nz+1)=q11(i,j,nz+1)*pvz(i,j,nz+1)
                qq2=pvx(i,j,nz+1)*pvx(i,j,nz+1)+pvy(i,j,nz+1)*pvy(i,j,nz+1)+pvz(i,j,nz+1)*pvz(i,j,nz+1)
                q15(i,j,nz+1)=2.5d0*p(i,j,nz+1)+0.5d0*q11(i,j,nz+1)*qq2
                q16(i,j,nz+1)=-q16(i,j,nz)
            end do
        end do
    end if
    !********处理通信数据，推出对应守恒量
    !*********轴向x向动静叶交界面
    if(slidm/=0)then
        if(rm==1)then
            i=nx+1
        else if(rm==2)then
            i=0
        end if
        q11(i,:,:)=0.
        pvx(i,:,:)=0.
        pvy(i,:,:)=0.
        pvz(i,:,:)=0.
        p(i,:,:)=0.
        q16(i,:,:)=0.
        do j=1,ny
            do k=1,nz
                do l=1,ssum
                    do ii=1,2
                        do jj=1,ny
                            do kk=bk(l),ek(l)
                                if(v(j,k,ii,jj,kk,l)/=0.)then
                                    q11(i,j,k)=q11(i,j,k)+v(j,k,ii,jj,kk,l)*sq11(ii,jj,kk,l)
                                    pvx(i,j,k)=pvx(i,j,k)+v(j,k,ii,jj,kk,l)*spvx(ii,jj,kk,l)
                                    pvy(i,j,k)=pvy(i,j,k)+v(j,k,ii,jj,kk,l)*spvy(ii,jj,kk,l)
                                    pvz(i,j,k)=pvz(i,j,k)+v(j,k,ii,jj,kk,l)*spvz(ii,jj,kk,l)
                                    p(i,j,k)=p(i,j,k)+v(j,k,ii,jj,kk,l)*sp(ii,jj,kk,l)
                                    q16(i,j,k)=q16(i,j,k)+v(j,k,ii,jj,kk,l)*sq16(ii,jj,kk,l)
                                end if
                            end do
                        end do
                    end do
                end do
                q12(i,j,k)= q11(i,j,k)*pvx(i,j,k)
                q13(i,j,k)= q11(i,j,k)*pvy(i,j,k)
                q14(i,j,k)= q11(i,j,k)*pvz(i,j,k)
                qq2=pvx(i,j,k)*pvx(i,j,k)+pvy(i,j,k)*pvy(i,j,k)+pvz(i,j,k)*pvz(i,j,k)
                q15(i,j,k)=2.5d0*p(i,j,k)+0.5d0*q11(i,j,k)*qq2
                t(i,j,k)=p(i,j,k)/(Q11(i,j,k)*rg)
            end do
        end do
        deallocate(sq11)
        deallocate(spvx)
        deallocate(spvy)
        deallocate(spvz)
        deallocate(sp)
        deallocate(sq16)
    end if
    !*********叶片内部轴向x向交界面
    if(xll>0)then
        i=0
        do k=1,nz
            do j=1,ny
                q12(i,j,k)= q11(i,j,k)*pvx(i,j,k)
                q13(i,j,k)= q11(i,j,k)*pvy(i,j,k)
                q14(i,j,k)= q11(i,j,k)*pvz(i,j,k)
                qq2=pvx(i,j,k)*pvx(i,j,k)+pvy(i,j,k)*pvy(i,j,k)+pvz(i,j,k)*pvz(i,j,k)
                q15(i,j,k)=2.5d0*p(i,j,k)+0.5d0*q11(i,j,k)*qq2
                t(i,j,k)=p(i,j,k)/(Q11(i,j,k)*rg)
            end do
        end do
    end if
    if(xll<xln(0)+xln(1)+xln(2)-1)then
        i=nx+1
        do k=1,nz
            do j=1,ny
                q12(i,j,k)= q11(i,j,k)*pvx(i,j,k)
                q13(i,j,k)= q11(i,j,k)*pvy(i,j,k)
                q14(i,j,k)= q11(i,j,k)*pvz(i,j,k)
                qq2=pvx(i,j,k)*pvx(i,j,k)+pvy(i,j,k)*pvy(i,j,k)+pvz(i,j,k)*pvz(i,j,k)
                q15(i,j,k)=2.5d0*p(i,j,k)+0.5d0*q11(i,j,k)*qq2
                t(i,j,k)=p(i,j,k)/(Q11(i,j,k)*rg)
            end do
        end do
    end if
    !*********y向叶片内交界面
    if(yl>0)then
        j=0
        do k=1,nz
            do i=1,nx
                q12(i,j,k)= q11(i,j,k)*pvx(i,j,k)
                q13(i,j,k)= q11(i,j,k)*pvy(i,j,k)
                q14(i,j,k)= q11(i,j,k)*pvz(i,j,k)
                qq2=pvx(i,j,k)*pvx(i,j,k)+pvy(i,j,k)*pvy(i,j,k)+pvz(i,j,k)*pvz(i,j,k)
                q15(i,j,k)=2.5d0*p(i,j,k)+0.5d0*q11(i,j,k)*qq2
                t(i,j,k)=p(i,j,k)/(Q11(i,j,k)*rg)
            end do
        end do
    end if
    if(yl<yln-1)then
        j=ny+1
        do k=1,nz
            do i=1,nx
                q12(i,j,k)= q11(i,j,k)*pvx(i,j,k)
                q13(i,j,k)= q11(i,j,k)*pvy(i,j,k)
                q14(i,j,k)= q11(i,j,k)*pvz(i,j,k)
                qq2=pvx(i,j,k)*pvx(i,j,k)+pvy(i,j,k)*pvy(i,j,k)+pvz(i,j,k)*pvz(i,j,k)
                q15(i,j,k)=2.5d0*p(i,j,k)+0.5d0*q11(i,j,k)*qq2
                t(i,j,k)=p(i,j,k)/(Q11(i,j,k)*rg)
            end do
        end do
    end if
    !*********周向z
    if(cycle_by<=ny)then
        do j=cycle_by,ny
            do i=1,nx
                q12(i,j,0)= q11(i,j,0)*pvx(i,j,0)
                q13(i,j,0)= q11(i,j,0)*pvy(i,j,0)
                q14(i,j,0)= q11(i,j,0)*pvz(i,j,0)
                qq2=pvx(i,j,0)*pvx(i,j,0)+pvy(i,j,0)*pvy(i,j,0)+pvz(i,j,0)*pvz(i,j,0)
                q15(i,j,0)=2.5d0*p(i,j,0)+0.5d0*q11(i,j,0)*qq2
                t(i,j,0)=p(i,j,0)/(Q11(i,j,0)*rg)
            
                q12(i,j,nz+1)= q11(i,j,nz+1)*pvx(i,j,nz+1)
                q13(i,j,nz+1)= q11(i,j,nz+1)*pvy(i,j,nz+1)
                q14(i,j,nz+1)= q11(i,j,nz+1)*pvz(i,j,nz+1)
                qq2=pvx(i,j,nz+1)*pvx(i,j,nz+1)+pvy(i,j,nz+1)*pvy(i,j,nz+1)+pvz(i,j,nz+1)*pvz(i,j,nz+1)
                q15(i,j,nz+1)=2.5d0*p(i,j,nz+1)+0.5d0*q11(i,j,nz+1)*qq2
                t(i,j,nz+1)=p(i,j,nz+1)/(Q11(i,j,nz+1)*rg)
            end do
        end do
    end if
    end subroutine bc
    
subroutine step !vol/t    !用佘程序（对流和粘性，谱半径也各向异性）
    use global
    implicit none
    real(8) :: tc,td,aa,cv,kc
    real(8),allocatable ::sri1(:,:,:),srj1(:,:,:),srk1(:,:,:) !修正后的谱半径（上加一横）
    
    do k=1,nz
        do j=1,ny
            do i=1,nx
                !******x方向
                vx=0.5d0*(pvx(i,j,k)+pvx(i+1,j,k))
                vy=0.5d0*(pvy(i,j,k)+pvy(i+1,j,k))
                vz=0.5d0*(pvz(i,j,k)+pvz(i+1,j,k))
                wx=vx
                wy=vy+rpm(rm)*zz02(i+1,j,k)
                wz=vz-rpm(rm)*yy02(i+1,j,k)
                tc=abs(wx*s2x(i+1,j,k)+wy*s2y(i+1,j,k)+wz*s2z(i+1,j,k))
                !******y方向
                vx=0.5d0*(pvx(i,j,k)+pvx(i,j+1,k))
                vy=0.5d0*(pvy(i,j,k)+pvy(i,j+1,k))
                vz=0.5d0*(pvz(i,j,k)+pvz(i,j+1,k))
                wx=vx
                wy=vy+rpm(rm)*zz03(i,j+1,k)
                wz=vz-rpm(rm)*yy03(i,j+1,k)
                tc=abs(wx*s3x(i,j+1,k)+wy*s3y(i,j+1,k)+wz*s3z(i,j+1,k))+tc
                !******z方向
                vx=0.5d0*(pvx(i,j,k)+pvx(i,j,k+1))
                vy=0.5d0*(pvy(i,j,k)+pvy(i,j,k+1))
                vz=0.5d0*(pvz(i,j,k)+pvz(i,j,k+1))
                wx=vx
                wy=vy+rpm(rm)*zz01(i,j,k+1)
                wz=vz-rpm(rm)*yy01(i,j,k+1)
                tc=abs(wx*s1x(i,j,k+1)+wy*s1y(i,j,k+1)+wz*s1z(i,j,k+1))+tc

                aa=sqrt(1.4d0*p(i,j,k)/q11(i,j,k))    !aa:音速c
                tc=tc+aa*(sqrt(s1x(i,j,k+1)**2+s1y(i,j,k+1)**2+s1z(i,j,k+1)**2)+sqrt(s2x(i+1,j,k)**2+&
                   s2y(i+1,j,k)**2+s2z(i+1,j,k)**2)+sqrt(s3x(i,j+1,k)**2+s3y(i,j+1,k)**2+s3z(i,j+1,k)**2))
                
                td=s1x(i,j,k+1)**2+s1y(i,j,k+1)**2+s1z(i,j,k+1)**2+s2x(i+1,j,k)**2+s2y(i+1,j,k)**2+&
                   s2z(i+1,j,k)**2+s3x(i,j+1,k)**2+s3y(i,j+1,k)**2+s3z(i,j+1,k)**2+2.d0*(abs(s1x(i,j,k+1)&
                   *s2x(i+1,j,k)+s1y(i,j,k+1)*s2y(i+1,j,k)+s1z(i,j,k+1)*s2z(i+1,j,k))+abs(s1x(i,j,k+1)&
                   *s3x(i,j+1,k)+s1y(i,j,k+1)*s3y(i,j+1,k)+s1z(i,j,k+1)*s3z(i,j+1,k))+abs(s3x(i,j+1,k)&
                   *s2x(i+1,j,k)+s3y(i,j+1,k)*s2y(i+1,j,k)+s3z(i,j+1,k)*s2z(i+1,j,k)))
                call viscosity(t(i,j,k),q16(i,j,k),cv,kc)
                td=8.d0*cv*td/q11(i,j,k)/vv0(i,j,k)
                time(i,j,k)=cfl*vv0(i,j,k)/(tc+td)
                sri(i,j,k)=tc+td
                srj(i,j,k)=tc+td
                srk(i,j,k)=tc+td
            end do
        end do
    end do
    !---------(needed by "ddd")各向异性lamt修正值，上加一横
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
    !——lamt修正值 赋值虚拟网格
    sri1(0,1:ny,1:nz) = sri1(1,1:ny,1:nz)
    srj1(0,1:ny,1:nz) = srj1(1,1:ny,1:nz)
    srk1(0,1:ny,1:nz) = srk1(1,1:ny,1:nz)
    sri1(nx+1,1:ny,1:nz) = sri1(nx,1:ny,1:nz)
    srj1(nx+1,1:ny,1:nz) = srj1(nx,1:ny,1:nz)
    srk1(nx+1,1:ny,1:nz) = srk1(nx,1:ny,1:nz)
    
    sri1(1:nx,0,1:nz) = sri1(1:nx,1,1:nz)
    srj1(1:nx,0,1:nz) = srj1(1:nx,1,1:nz)
    srk1(1:nx,0,1:nz) = srk1(1:nx,1,1:nz)
    sri1(1:nx,ny+1,1:nz) = sri1(1:nx,ny,1:nz)
    srj1(1:nx,ny+1,1:nz) = srj1(1:nx,ny,1:nz)
    srk1(1:nx,ny+1,1:nz) = srk1(1:nx,ny,1:nz)
            
    sri1(1:nx,1:ny,0) = sri1(1:nx,1:ny,1)
    srj1(1:nx,1:ny,0) = srj1(1:nx,1:ny,1)
    srk1(1:nx,1:ny,0) = srk1(1:nx,1:ny,1)
    sri1(1:nx,1:ny,nz+1) = sri1(1:nx,1:ny,nz)
    srj1(1:nx,1:ny,nz+1) = srj1(1:nx,1:ny,nz)
    srk1(1:nx,1:ny,nz+1) = srk1(1:nx,1:ny,nz)
    !——lamt修正值 赋值lamt
    sri(0:nx+1,0:ny+1,0:nz+1) = sri1(0:nx+1,0:ny+1,0:nz+1)
    srj(0:nx+1,0:ny+1,0:nz+1) = srj1(0:nx+1,0:ny+1,0:nz+1)
    srk(0:nx+1,0:ny+1,0:nz+1) = srk1(0:nx+1,0:ny+1,0:nz+1)
    deallocate(sri1)
    deallocate(srj1)
    deallocate(srk1)
    end subroutine step
    
subroutine ddd !需要先调用当地时间步长             !xyz全部包含边界
    use global
    implicit none
    real(8) ::em2,em4 !二阶、四阶系数
    real(8),allocatable ::dp(:)     !压力二阶导，两个方向上，定义一个动态数组，两个分别前后用
    real(8) ::flu1,flu2,flu3,flu4,flu5,flu6,ram
    
    av1=0.d0
    av2=0.d0
    av3=0.d0
    av4=0.d0
    av5=0.d0
    av6=0.d0
    q15(0:nx+1,0:ny+1,0:nz+1)=q15(0:nx+1,0:ny+1,0:nz+1)+p(0:nx+1,0:ny+1,0:nz+1)  !rou*h=rou*e+p
    ! i-direction -----------------------------------------------------------------
    allocate(dp(0:nx+1))
    dp(0)=0.d0
    dp(nx+1)=0.d0
    do k=1,nz  !沿x方向的压力二阶导
      do j=1,ny
         do i=1,nx
             dp(i)=abs((p(i+1,j,k)-2.d0*p(i,j,k)+p(i-1,j,k))/(p(i+1,j,k)+2.d0*p(i,j,k)+p(i-1,j,k)))
         end do
     ! - dissipation fluxes (at I+1/2)
        do i=0,nx
           ram   = 0.5d0*(sri(i,j,k)+sri(i+1,j,k))!lamt
           em2   =a2*ram*max(dp(i),dp(i+1))
           em4   =a4*ram
           em4   =max(em4-em2,0.d0)       !求em4-em2和0的最大值，这里的em2,em4都是公式前差算子的总系数。
           
           flu1  =em2*(Q11(i+1,j,k)-Q11(i,j,k))
           flu2  =em2*(Q12(i+1,j,k)-Q12(i,j,k))
           flu3  =em2*(Q13(i+1,j,k)-Q13(i,j,k))
           flu4  =em2*(Q14(i+1,j,k)-Q14(i,j,k))
           flu5  =em2*(Q15(i+1,j,k)-Q15(i,j,k))
           flu6  =em2*(Q16(i+1,j,k)-Q16(i,j,k))
           
           if(i>=1 .and. i<=nx-1) then        !改成1，nx-1后算到规定步数还没收敛到-3，算了6000步才收敛到-2.1
               flu1=flu1+em4*(q11(i-1,j,k)-3.d0*q11(i,j,k)+3.d0*q11(i+1,j,k)-q11(i+2,j,k))
               flu2=flu2+em4*(q12(i-1,j,k)-3.d0*q12(i,j,k)+3.d0*q12(i+1,j,k)-q12(i+2,j,k))
               flu3=flu3+em4*(q13(i-1,j,k)-3.d0*q13(i,j,k)+3.d0*q13(i+1,j,k)-q13(i+2,j,k))
               flu4=flu4+em4*(q14(i-1,j,k)-3.d0*q14(i,j,k)+3.d0*q14(i+1,j,k)-q14(i+2,j,k))
               flu5=flu5+em4*(q15(i-1,j,k)-3.d0*q15(i,j,k)+3.d0*q15(i+1,j,k)-q15(i+2,j,k))
               flu6=flu6+em4*(q16(i-1,j,k)-3.d0*q16(i,j,k)+3.d0*q16(i+1,j,k)-q16(i+2,j,k))
           end if
           ! --- dissipation term
           av1(i  ,j,k) = av1(i  ,j,k) +flu1   !人工粘性通量flu在网格的右侧面，右侧流出为正，左侧流入为负。
           av2(i  ,j,k) = av2(i  ,j,k) +flu2
           av3(i  ,j,k) = av3(i  ,j,k) +flu3
           av4(i  ,j,k) = av4(i  ,j,k) +flu4
           av5(i  ,j,k) = av5(i  ,j,k) +flu5
           av6(i  ,j,k) = av6(i  ,j,k) +flu6
          
           av1(i+1,j,k) = av1(i+1,j,k) -flu1
           av2(i+1,j,k) = av2(i+1,j,k) -flu2
           av3(i+1,j,k) = av3(i+1,j,k) -flu3
           av4(i+1,j,k) = av4(i+1,j,k) -flu4
           av5(i+1,j,k) = av5(i+1,j,k) -flu5
           av6(i+1,j,k) = av6(i+1,j,k) -flu6
        end do
      end do
    end do
    deallocate(dp)
    ! j-direction -----------------------------------------------------------------
    allocate(dp(0:ny+1))
    dp(0)=0.d0
    dp(ny+1)=0.d0
    do k=1,nz  !沿x方向的压力二阶导
      do i=1,nx
         do j=1,ny
             dp(j)=abs((p(i,j+1,k)-2.d0*p(i,j,k)+p(i,j-1,k))/(p(i,j+1,k)+2.d0*p(i,j,k)+p(i,j-1,k)))
         end do
     ! - dissipation fluxes (at I+1/2)
        do j=0,ny
           ram   = 0.5d0*(srj(i,j,k)+srj(i,j+1,k))!lamt
           em2   =a2*ram*max(dp(j),dp(j+1))
           em4   =a4*ram
           em4   =max(em4-em2,0.d0)       !求em4-em2和0的最大值，这里的em2,em4都是公式前差算子的总系数。
           
           flu1  =em2*(Q11(i,j+1,k)-Q11(i,j,k))
           flu2  =em2*(Q12(i,j+1,k)-Q12(i,j,k))
           flu3  =em2*(Q13(i,j+1,k)-Q13(i,j,k))
           flu4  =em2*(Q14(i,j+1,k)-Q14(i,j,k))
           flu5  =em2*(Q15(i,j+1,k)-Q15(i,j,k))
           flu6  =em2*(Q16(i,j+1,k)-Q16(i,j,k))
            
           if(j>=1 .and. j<=ny-1) then        !改成1，nx-1后算到规定步数还没收敛到-3，算了6000步才收敛到-2.1
               flu1=flu1+em4*(q11(i,j-1,k)-3.d0*q11(i,j,k)+3.d0*q11(i,j+1,k)-q11(i,j+2,k))
               flu2=flu2+em4*(q12(i,j-1,k)-3.d0*q12(i,j,k)+3.d0*q12(i,j+1,k)-q12(i,j+2,k))
               flu3=flu3+em4*(q13(i,j-1,k)-3.d0*q13(i,j,k)+3.d0*q13(i,j+1,k)-q13(i,j+2,k))
               flu4=flu4+em4*(q14(i,j-1,k)-3.d0*q14(i,j,k)+3.d0*q14(i,j+1,k)-q14(i,j+2,k))
               flu5=flu5+em4*(q15(i,j-1,k)-3.d0*q15(i,j,k)+3.d0*q15(i,j+1,k)-q15(i,j+2,k))
               flu6=flu6+em4*(q16(i,j-1,k)-3.d0*q16(i,j,k)+3.d0*q16(i,j+1,k)-q16(i,j+2,k))
           end if
           ! --- dissipation term
           av1(i  ,j,k) = av1(i  ,j,k) +flu1   !人工粘性通量flu在网格的右侧面，右侧流出为正，左侧流入为负。
           av2(i  ,j,k) = av2(i  ,j,k) +flu2
           av3(i  ,j,k) = av3(i  ,j,k) +flu3
           av4(i  ,j,k) = av4(i  ,j,k) +flu4
           av5(i  ,j,k) = av5(i  ,j,k) +flu5
           av6(i  ,j,k) = av6(i  ,j,k) +flu6
           
           av1(i,j+1,k) = av1(i,j+1,k) -flu1
           av2(i,j+1,k) = av2(i,j+1,k) -flu2
           av3(i,j+1,k) = av3(i,j+1,k) -flu3
           av4(i,j+1,k) = av4(i,j+1,k) -flu4
           av5(i,j+1,k) = av5(i,j+1,k) -flu5
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
                  ram   = 0.5d0*(srk(i,j,k)+srk(i,j,k+1))!lamt
                  em2   =a2*ram*max(dp(k),dp(k+1))
                  em4   =a4*ram
                  em4   =max(em4-em2,0.d0)       !求em4-em2和0的最大值，这里的em2,em4都是公式前差算子的总系数。
           
                  flu1  =em2*(Q11(i,j,k+1)-Q11(i,j,k))
                  flu2  =em2*(Q12(i,j,k+1)-Q12(i,j,k))
                  flu3  =em2*(Q13(i,j,k+1)-Q13(i,j,k))
                  flu4  =em2*(Q14(i,j,k+1)-Q14(i,j,k))
                  flu5  =em2*(Q15(i,j,k+1)-Q15(i,j,k))
                  flu6  =em2*(Q16(i,j,k+1)-Q16(i,j,k))
                  
                  if(k>=1 .and. k<=nz-1) then        !改成1，nx-1后算到规定步数还没收敛到-3，算了6000步才收敛到-2.1
                      flu1=flu1+em4*(q11(i,j,k-1)-3.d0*q11(i,j,k)+3.d0*q11(i,j,k+1)-q11(i,j,k+2))
                      flu2=flu2+em4*(q12(i,j,k-1)-3.d0*q12(i,j,k)+3.d0*q12(i,j,k+1)-q12(i,j,k+2))
                      flu3=flu3+em4*(q13(i,j,k-1)-3.d0*q13(i,j,k)+3.d0*q13(i,j,k+1)-q13(i,j,k+2))
                      flu4=flu4+em4*(q14(i,j,k-1)-3.d0*q14(i,j,k)+3.d0*q14(i,j,k+1)-q14(i,j,k+2))
                      flu5=flu5+em4*(q15(i,j,k-1)-3.d0*q15(i,j,k)+3.d0*q15(i,j,k+1)-q15(i,j,k+2))
                      flu6=flu6+em4*(q16(i,j,k-1)-3.d0*q16(i,j,k)+3.d0*q16(i,j,k+1)-q16(i,j,k+2))
                  end if
                  ! --- dissipation term
                  av1(i  ,j,k) = av1(i  ,j,k) +flu1   !人工粘性通量flu在网格的右侧面，右侧流出为正，左侧流入为负。
                  av2(i  ,j,k) = av2(i  ,j,k) +flu2
                  av3(i  ,j,k) = av3(i  ,j,k) +flu3
                  av4(i  ,j,k) = av4(i  ,j,k) +flu4
                  av5(i  ,j,k) = av5(i  ,j,k) +flu5
                  av6(i  ,j,k) = av6(i  ,j,k) +flu6
                  
                  av1(i,j,k+1) = av1(i,j,k+1) -flu1
                  av2(i,j,k+1) = av2(i,j,k+1) -flu2
                  av3(i,j,k+1) = av3(i,j,k+1) -flu3
                  av4(i,j,k+1) = av4(i,j,k+1) -flu4
                  av5(i,j,k+1) = av5(i,j,k+1) -flu5
                  av6(i,j,k+1) = av6(i,j,k+1) -flu6
               end do
       end do
    end do
    q15(0:nx+1,0:ny+1,0:nz+1)=q15(0:nx+1,0:ny+1,0:nz+1)-p(0:nx+1,0:ny+1,0:nz+1)
    end subroutine ddd

subroutine qqq
    use global
    implicit none
    real(8) ::flu1,flu2,flu3,flu4,flu5,flu6,tur,vf,rf !界面上平均值
    
    qc1=0.d0
    qc2=0.d0
    qc3=0.d0
    qc4=0.d0
    qc5=0.d0
    qc6=0.d0
    !******************x方向********************
    do k=1,nz
        do j=1,ny
            do i=1,nx+1
               vx=0.5d0*(pvx(i,j,k)+pvx(i-1,j,k))
               vy=0.5d0*(pvy(i,j,k)+pvy(i-1,j,k))
               vz=0.5d0*(pvz(i,j,k)+pvz(i-1,j,k))
               vf=vx*s2x(i,j,k)+vy*s2y(i,j,k)+vz*s2z(i,j,k)             !绝对速度通量
               rf=rpm(rm)*zz02(i,j,k)*s2y(i,j,k)-rpm(rm)*yy02(i,j,k)*s2z(i,j,k)     !牵连速度通量
               vf=vf+rf                                                 !相对速度通量
               dim=0.5d0*(Q11(i,j,k)+Q11(i-1,j,k))
               pp=0.5d0*(p(i,j,k)+p(i-1,j,k))
               en=0.5d0*(Q15(i,j,k)+Q15(i-1,j,k))
               tur=0.5d0*(Q16(i,j,k)+Q16(i-1,j,k))
               
               flu1=dim*vf
               flu2=flu1*vx+pp*s2x(i,j,k)
               flu3=flu1*vy+pp*s2y(i,j,k)
               flu4=flu1*vz+pp*s2z(i,j,k)
               flu5=(en+pp)*vf-pp*rf
               flu6=tur*vf
               
               qc1(i,j,k)=qc1(i,j,k)+flu1
               qc2(i,j,k)=qc2(i,j,k)+flu2
               qc3(i,j,k)=qc3(i,j,k)+flu3
               qc4(i,j,k)=qc4(i,j,k)+flu4
               qc5(i,j,k)=qc5(i,j,k)+flu5
               qc6(i,j,k)=qc6(i,j,k)+flu6
               
               qc1(i-1,j,k)=qc1(i-1,j,k)-flu1
               qc2(i-1,j,k)=qc2(i-1,j,k)-flu2
               qc3(i-1,j,k)=qc3(i-1,j,k)-flu3
               qc4(i-1,j,k)=qc4(i-1,j,k)-flu4
               qc5(i-1,j,k)=qc5(i-1,j,k)-flu5
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
               rf=rpm(rm)*zz03(i,j,k)*s3y(i,j,k)-rpm(rm)*yy03(i,j,k)*s3z(i,j,k)     !牵连速度通量
               vf=vf+rf                                                 !相对速度通量
               dim=0.5d0*(Q11(i,j,k)+Q11(i,j-1,k))
               pp=0.5d0*(p(i,j,k)+p(i,j-1,k))
               en=0.5d0*(Q15(i,j,k)+Q15(i,j-1,k))
               tur=0.5d0*(Q16(i,j,k)+Q16(i,j-1,k))
               
               flu1=dim*vf
               flu2=flu1*vx+pp*s3x(i,j,k)
               flu3=flu1*vy+pp*s3y(i,j,k)
               flu4=flu1*vz+pp*s3z(i,j,k)
               flu5=(en+pp)*vf-pp*rf
               flu6=tur*vf
               
               qc1(i,j,k)=qc1(i,j,k)+flu1
               qc2(i,j,k)=qc2(i,j,k)+flu2
               qc3(i,j,k)=qc3(i,j,k)+flu3
               qc4(i,j,k)=qc4(i,j,k)+flu4
               qc5(i,j,k)=qc5(i,j,k)+flu5
               qc6(i,j,k)=qc6(i,j,k)+flu6
               
               qc1(i,j-1,k)=qc1(i,j-1,k)-flu1
               qc2(i,j-1,k)=qc2(i,j-1,k)-flu2
               qc3(i,j-1,k)=qc3(i,j-1,k)-flu3
               qc4(i,j-1,k)=qc4(i,j-1,k)-flu4
               qc5(i,j-1,k)=qc5(i,j-1,k)-flu5
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
                rf=rpm(rm)*zz01(i,j,k)*s1y(i,j,k)-rpm(rm)*yy01(i,j,k)*s1z(i,j,k)     !牵连速度通量
                vf=vf+rf                                                 !相对速度通量
                dim=0.5d0*(Q11(i,j,k)+Q11(i,j,k-1))
                pp=0.5d0*(p(i,j,k)+p(i,j,k-1))
                en=0.5d0*(Q15(i,j,k)+Q15(i,j,k-1))
                tur=0.5d0*(Q16(i,j,k)+Q16(i,j,k-1))
                
                flu1=dim*vf
                flu2=flu1*vx+pp*s1x(i,j,k)
                flu3=flu1*vy+pp*s1y(i,j,k)
                flu4=flu1*vz+pp*s1z(i,j,k)
                flu5=(en+pp)*vf-pp*rf
                flu6=tur*vf
                
                qc1(i,j,k)=qc1(i,j,k)+flu1
                qc2(i,j,k)=qc2(i,j,k)+flu2
                qc3(i,j,k)=qc3(i,j,k)+flu3
                qc4(i,j,k)=qc4(i,j,k)+flu4
                qc5(i,j,k)=qc5(i,j,k)+flu5
                qc6(i,j,k)=qc6(i,j,k)+flu6
                
                qc1(i,j,k-1)=qc1(i,j,k-1)-flu1
                qc2(i,j,k-1)=qc2(i,j,k-1)-flu2
                qc3(i,j,k-1)=qc3(i,j,k-1)-flu3
                qc4(i,j,k-1)=qc4(i,j,k-1)-flu4
                qc5(i,j,k-1)=qc5(i,j,k-1)-flu5
                qc6(i,j,k-1)=qc6(i,j,k-1)-flu6
            end do
        end do
    end do
    end subroutine qqq
    
subroutine qqqv
    use global
    implicit none
    real(8) ::flu2,flu3,flu4,flu5,flu6,tauxx,tauyy,tauzz,tauxy,tauxz,tauyz,phix,phiy,phiz
    real(8) ::uav,vav,wav,q16av,tav,mav,kav
    real(8),allocatable ::tur(:,:,:),gradfi(:,:,:,:),gradfj(:,:,:,:),gradfk(:,:,:,:)
    
    allocate(tur(0:nx+1,0:ny+1,0:nz+1))
    tur(0:nx+1,0:ny+1,0:nz+1)=q16(0:nx+1,0:ny+1,0:nz+1)/q11(0:nx+1,0:ny+1,0:nz+1)
    !计算界面上的导数
    allocate(gradfi(15,1:nx+1,0:ny+1,0:nz+1))
    allocate(gradfj(15,0:nx+1,1:ny+1,0:nz+1))
    allocate(gradfk(15,0:nx+1,0:ny+1,1:nz+1))
    call gradsfaceI(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfi(1,1:nx+1,0:ny+1,0:nz+1)) !du/dx
    call gradsfaceI(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfi(2,1:nx+1,0:ny+1,0:nz+1)) !du/dy
    call gradsfaceI(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfi(3,1:nx+1,0:ny+1,0:nz+1)) !du/dz
    call gradsfaceI(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfi(4,1:nx+1,0:ny+1,0:nz+1)) !dv/dx
    call gradsfaceI(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfi(5,1:nx+1,0:ny+1,0:nz+1)) !dv/dy
    call gradsfaceI(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfi(6,1:nx+1,0:ny+1,0:nz+1)) !dv/dz
    call gradsfaceI(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfi(7,1:nx+1,0:ny+1,0:nz+1)) !dw/dx
    call gradsfaceI(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfi(8,1:nx+1,0:ny+1,0:nz+1)) !dw/dy
    call gradsfaceI(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfi(9,1:nx+1,0:ny+1,0:nz+1)) !dw/dz
    call gradsfaceI(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfi(10,1:nx+1,0:ny+1,0:nz+1)) !dt/dx
    call gradsfaceI(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfi(11,1:nx+1,0:ny+1,0:nz+1)) !dt/dy
    call gradsfaceI(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfi(12,1:nx+1,0:ny+1,0:nz+1)) !dt/dz
    call gradsfaceI(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfi(13,1:nx+1,0:ny+1,0:nz+1)) !dtur/dx
    call gradsfaceI(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfi(14,1:nx+1,0:ny+1,0:nz+1)) !dtur/dy
    call gradsfaceI(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfi(15,1:nx+1,0:ny+1,0:nz+1)) !dtur/dz

    call gradsfaceJ(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfj(1,0:nx+1,1:ny+1,0:nz+1)) !du/dx
    call gradsfaceJ(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfj(2,0:nx+1,1:ny+1,0:nz+1)) !du/dy
    call gradsfaceJ(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfj(3,0:nx+1,1:ny+1,0:nz+1)) !du/dz
    call gradsfaceJ(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfj(4,0:nx+1,1:ny+1,0:nz+1)) !dv/dx
    call gradsfaceJ(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfj(5,0:nx+1,1:ny+1,0:nz+1)) !dv/dy
    call gradsfaceJ(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfj(6,0:nx+1,1:ny+1,0:nz+1)) !dv/dz
    call gradsfaceJ(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfj(7,0:nx+1,1:ny+1,0:nz+1)) !dw/dx
    call gradsfaceJ(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfj(8,0:nx+1,1:ny+1,0:nz+1)) !dw/dy
    call gradsfaceJ(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfj(9,0:nx+1,1:ny+1,0:nz+1)) !dw/dz
    call gradsfaceJ(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfj(10,0:nx+1,1:ny+1,0:nz+1)) !dt/dx
    call gradsfaceJ(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfj(11,0:nx+1,1:ny+1,0:nz+1)) !dt/dy
    call gradsfaceJ(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfj(12,0:nx+1,1:ny+1,0:nz+1)) !dt/dz
    call gradsfaceJ(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfj(13,0:nx+1,1:ny+1,0:nz+1)) !dtur/dx
    call gradsfaceJ(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfj(14,0:nx+1,1:ny+1,0:nz+1)) !dtur/dy
    call gradsfaceJ(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfj(15,0:nx+1,1:ny+1,0:nz+1)) !dtur/dz

    call gradsfaceK(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfk(1,0:nx+1,0:ny+1,1:nz+1)) !du/dx
    call gradsfaceK(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfk(2,0:nx+1,0:ny+1,1:nz+1)) !du/dy
    call gradsfaceK(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfk(3,0:nx+1,0:ny+1,1:nz+1)) !du/dz
    call gradsfaceK(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfk(4,0:nx+1,0:ny+1,1:nz+1)) !dv/dx
    call gradsfaceK(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfk(5,0:nx+1,0:ny+1,1:nz+1)) !dv/dy
    call gradsfaceK(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvy(0:nx+1,0:ny+1,0:nz+1),gradfk(6,0:nx+1,0:ny+1,1:nz+1)) !dv/dz
    call gradsfaceK(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfk(7,0:nx+1,0:ny+1,1:nz+1)) !dw/dx
    call gradsfaceK(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfk(8,0:nx+1,0:ny+1,1:nz+1)) !dw/dy
    call gradsfaceK(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),pvz(0:nx+1,0:ny+1,0:nz+1),gradfk(9,0:nx+1,0:ny+1,1:nz+1)) !dw/dz
    call gradsfaceK(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfk(10,0:nx+1,0:ny+1,1:nz+1)) !dt/dx
    call gradsfaceK(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfk(11,0:nx+1,0:ny+1,1:nz+1)) !dt/dy
    call gradsfaceK(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),t(0:nx+1,0:ny+1,0:nz+1),gradfk(12,0:nx+1,0:ny+1,1:nz+1)) !dt/dz
    call gradsfaceK(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfk(13,0:nx+1,0:ny+1,1:nz+1)) !dtur/dx
    call gradsfaceK(s2y(0:nx+2,1:ny,1:nz),s3y(1:nx,0:ny+2,1:nz),s1y(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfk(14,0:nx+1,0:ny+1,1:nz+1)) !dtur/dy
    call gradsfaceK(s2z(0:nx+2,1:ny,1:nz),s3z(1:nx,0:ny+2,1:nz),s1z(1:nx,1:ny,0:nz+2),tur(0:nx+1,0:ny+1,0:nz+1),gradfk(15,0:nx+1,0:ny+1,1:nz+1)) !dtur/dz

    qv2=0.d0
    qv3=0.d0
    qv4=0.d0
    qv5=0.d0
    qv6=0.d0
    ! i-direction -----------------------------------------------------------------
    do k=1,nz
        do j=1,ny
            do i=1,nx+1
                uav   = 0.5D0*(pvx(i-1,j,k)+pvx(i,j,k))     !当前网格左侧面量值
                vav   = 0.5D0*(pvy(i-1,j,k)+pvy(i,j,k))
                wav   = 0.5D0*(pvz(i-1,j,k)+pvz(i,j,k))
                q16av = 0.5D0*(q16(i-1,j,k)+q16(i,j,k))
                tav   = 0.5D0*(t(i-1,j,k)+t(i,j,k))
                cvl   = cvl0*(tav/t0)**1.5*(t0+ts)/(tav+ts)
                call viscosity(tav,q16av,mav,kav)               !计算粘性系数和导热系数
                tauxx = two3*mav*(2.D0*gradfi(1,i,j,k)-gradfi(5,i,j,k)-gradfi(9,i,j,k))!沿x,z方向求时，对应的tauxx都是不一样的，要分别再求。
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
                flu6  = (cvl+q16av)*(gradfi(13,i,j,k)*s2x(i,j,k)+gradfi(14,i,j,k)*s2y(i,j,k)+gradfi(15,i,j,k)*s2z(i,j,k))/sigmav
            ! --- dissipation term
                qv2(i,j,k)=qv2(i,j,k)+flu2
                qv3(i,j,k)=qv3(i,j,k)+flu3
                qv4(i,j,k)=qv4(i,j,k)+flu4      !左侧流入为正，右侧流出为负；不用管粘性要负，总表达式前有负号
                qv5(i,j,k)=qv5(i,j,k)+flu5
                qv6(i,j,k) =qv6(i,j,k)+flu6

                qv2(i-1,j,k)=qv2(i-1,j,k)-flu2
                qv3(i-1,j,k)=qv3(i-1,j,k)-flu3
                qv4(i-1,j,k)=qv4(i-1,j,k)-flu4
                qv5(i-1,j,k)=qv5(i-1,j,k)-flu5
                qv6(i-1,j,k)=qv6(i-1,j,k)-flu6
            end do
        end do
    end do
    !j-direction -----------------------------------------------------------------
    do k=1,nz
        do j=1,ny+1
            do i=1,nx
                uav   = 0.5D0*(pvx(i,j-1,k)+pvx(i,j,k))     !当前网格左侧面量值
                vav   = 0.5D0*(pvy(i,j-1,k)+pvy(i,j,k))
                wav   = 0.5D0*(pvz(i,j-1,k)+pvz(i,j,k))
                q16av = 0.5D0*(q16(i,j-1,k)+q16(i,j,k))
                tav   = 0.5D0*(t(i,j-1,k)+t(i,j,k))
                cvl   = cvl0*(tav/t0)**1.5*(t0+ts)/(tav+ts)
                call viscosity(tav,q16av,mav,kav)               !计算粘性系数和导热系数
                tauxx = two3*mav*(2.D0*gradfj(1,i,j,k)-gradfj(5,i,j,k)-gradfj(9,i,j,k))!沿x,z方向求时，对应的tauxx都是不一样的，要分别再求。
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
                flu6  = (cvl+q16av)*(gradfj(13,i,j,k)*s3x(i,j,k)+gradfj(14,i,j,k)*s3y(i,j,k)+gradfj(15,i,j,k)*s3z(i,j,k))/sigmav
            ! --- dissipation term
                qv2(i,j,k)=qv2(i,j,k)+flu2
                qv3(i,j,k)=qv3(i,j,k)+flu3
                qv4(i,j,k)=qv4(i,j,k)+flu4      !左侧流入为正，右侧流出为负；不用管粘性要负，总表达式前有负号
                qv5(i,j,k)=qv5(i,j,k)+flu5
                qv6(i,j,k) =qv6(i,j,k)+flu6

                qv2(i,j-1,k)=qv2(i,j-1,k)-flu2
                qv3(i,j-1,k)=qv3(i,j-1,k)-flu3
                qv4(i,j-1,k)=qv4(i,j-1,k)-flu4
                qv5(i,j-1,k)=qv5(i,j-1,k)-flu5
                qv6(i,j-1,k)=qv6(i,j-1,k)-flu6
            end do
        end do
    end do
    ! k-direction -----------------------------------------------------------------
    do k=1,nz+1
        do j=1,ny
            do i=1,nx
                uav   = 0.5D0*(pvx(i,j,k-1)+pvx(i,j,k))     !当前网格左侧面量值
                vav   = 0.5D0*(pvy(i,j,k-1)+pvy(i,j,k))
                wav   = 0.5D0*(pvz(i,j,k-1)+pvz(i,j,k))
                q16av = 0.5D0*(q16(i,j,k-1)+q16(i,j,k))
                tav   = 0.5D0*(t(i,j,k-1)+t(i,j,k))
                cvl   = cvl0*(tav/t0)**1.5*(t0+ts)/(tav+ts)
                call viscosity(tav,q16av,mav,kav)               !计算粘性系数和导热系数
                tauxx = two3*mav*(2.D0*gradfk(1,i,j,k)-gradfk(5,i,j,k)-gradfk(9,i,j,k))!沿x,z方向求时，对应的tauxx都是不一样的，要分别再求。
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
                flu6  = (cvl+q16av)*(gradfk(13,i,j,k)*s1x(i,j,k)+gradfk(14,i,j,k)*s1y(i,j,k)+gradfk(15,i,j,k)*s1z(i,j,k))/sigmav
            ! --- dissipation term
                qv2(i,j,k)=qv2(i,j,k)+flu2
                qv3(i,j,k)=qv3(i,j,k)+flu3
                qv4(i,j,k)=qv4(i,j,k)+flu4      !左侧流入为正，右侧流出为负；不用管粘性要负，总表达式前有负号
                qv5(i,j,k)=qv5(i,j,k)+flu5
                qv6(i,j,k) =qv6(i,j,k)+flu6
                
                qv2(i,j,k-1)=qv2(i,j,k-1)-flu2
                qv3(i,j,k-1)=qv3(i,j,k-1)-flu3
                qv4(i,j,k-1)=qv4(i,j,k-1)-flu4
                qv5(i,j,k-1)=qv5(i,j,k-1)-flu5
                qv6(i,j,k-1)=qv6(i,j,k-1)-flu6
            end do
        end do
    end do
    deallocate(tur)
    deallocate(gradfi)
    deallocate(gradfj)
    deallocate(gradfk)
    !源项里的加到QV上 -----------------------------------------------------------------
    qv3(1:nx,1:ny,1:nz)=qv3(1:nx,1:ny,1:nz)+rpm(rm)*q14(1:nx,1:ny,1:nz)*vv0(1:nx,1:ny,1:nz)
    qv4(1:nx,1:ny,1:nz)=qv4(1:nx,1:ny,1:nz)-rpm(rm)*q13(1:nx,1:ny,1:nz)*vv0(1:nx,1:ny,1:nz)
    call SAsource
    end subroutine qqqv

subroutine SAsource
    use global
    implicit none
    real(8),allocatable ::tur(:,:,:),pwx(:,:,:),pwy(:,:,:),pwz(:,:,:),gradc(:,:,:,:),gradcs(:,:,:,:)
    real(8) ::tem,fv1,fv2,fv3,rp1,rp2,rp3,w12,w13,w23,w21,w31,w32,ww2,ww,svot,vm,gv
    real(8) ::s11,s22,s33,s12,s13,s23,ss2,sss,dd,fr1r1,fr1r2,fr1,yv,ra,ga,fw,dv

    allocate(tur(0:nx+1,0:ny+1,0:nz+1))
    allocate(pwx(0:nx+1,0:ny+1,0:nz+1))
    allocate(pwy(0:nx+1,0:ny+1,0:nz+1))
    allocate(pwz(0:nx+1,0:ny+1,0:nz+1))
    tur(0:nx+1,0:ny+1,0:nz+1)=q16(0:nx+1,0:ny+1,0:nz+1)/q11(0:nx+1,0:ny+1,0:nz+1)
    pwx(0:nx+1,0:ny+1,0:nz+1)=pvx(0:nx+1,0:ny+1,0:nz+1)
    pwy(0:nx+1,0:ny+1,0:nz+1)=pvy(0:nx+1,0:ny+1,0:nz+1)+rpm(rm)*zz0(0:nx+1,0:ny+1,0:nz+1)
    pwz(0:nx+1,0:ny+1,0:nz+1)=pvz(0:nx+1,0:ny+1,0:nz+1)-rpm(rm)*yy0(0:nx+1,0:ny+1,0:nz+1)
    !********网格中心导数***********
    allocate(gradc(12,0:nx+1,0:ny+1,0:nz+1))
    call gradscentre(1,pwx(0:nx+1,0:ny+1,0:nz+1),gradc(1,0:nx+1,0:ny+1,0:nz+1)) !du/dx
    call gradscentre(2,pwx(0:nx+1,0:ny+1,0:nz+1),gradc(2,0:nx+1,0:ny+1,0:nz+1)) !du/dy
    call gradscentre(3,pwx(0:nx+1,0:ny+1,0:nz+1),gradc(3,0:nx+1,0:ny+1,0:nz+1)) !du/dz
    call gradscentre(1,pwy(0:nx+1,0:ny+1,0:nz+1),gradc(4,0:nx+1,0:ny+1,0:nz+1)) !dw/dx
    call gradscentre(2,pwy(0:nx+1,0:ny+1,0:nz+1),gradc(5,0:nx+1,0:ny+1,0:nz+1)) !dw/dy
    call gradscentre(3,pwy(0:nx+1,0:ny+1,0:nz+1),gradc(6,0:nx+1,0:ny+1,0:nz+1)) !dw/dz
    call gradscentre(1,pwz(0:nx+1,0:ny+1,0:nz+1),gradc(7,0:nx+1,0:ny+1,0:nz+1)) !dw/dx
    call gradscentre(2,pwz(0:nx+1,0:ny+1,0:nz+1),gradc(8,0:nx+1,0:ny+1,0:nz+1)) !dw/dy
    call gradscentre(3,pwz(0:nx+1,0:ny+1,0:nz+1),gradc(9,0:nx+1,0:ny+1,0:nz+1)) !dw/dz
    call gradscentre(1,tur(0:nx+1,0:ny+1,0:nz+1),gradc(10,0:nx+1,0:ny+1,0:nz+1))  !dptur/dx
    call gradscentre(2,tur(0:nx+1,0:ny+1,0:nz+1),gradc(11,0:nx+1,0:ny+1,0:nz+1))  !dptur/dy
    call gradscentre(3,tur(0:nx+1,0:ny+1,0:nz+1),gradc(12,0:nx+1,0:ny+1,0:nz+1))  !dptur/dz
    !求r~里的ds/dt
    allocate(gradcs(9,0:nx+1,0:ny+1,0:nz+1))
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
                !********扩散项计算***********
                dv=cb2/sigmav*(gradc(10,i,j,k)**2+gradc(11,i,j,k)**2+gradc(12,i,j,k)**2)*q11(i,j,k)
                !********生成项gv计算***********
                cvl=cvl0*(t(i,j,k)/t0)**1.5*(t0+ts)/(t(i,j,k)+ts)
                tem=q16(i,j,k)/cvl
                tem=max(tem,1.d-4)
                fv1=1.d0/(1.d0+(cv1/tem)**3)
                fv2=(1.d0+tem/cv2)**(-3)
                fv3=(1.d0+tem*fv1)*(1.d0-fv2)/tem
                w12=0.5d0*(gradc(2,i,j,k)-gradc(4,i,j,k))
                w13=0.5d0*(gradc(3,i,j,k)-gradc(7,i,j,k))
                w23=0.5d0*(gradc(6,i,j,k)-gradc(8,i,j,k))
                ww2=4.d0*(w12*w12+w13*w13+w23*w23)
                ww=sqrt(ww2)               !论文里红圈圈
                vm=tur(i,j,k)/(kap*kap*dmini(i,j,k)*dmini(i,j,k))
                svot=ww*fv3+vm*fv2
                if (svot<0.d0) then
                    write(*,*) 'svot<0',svot,fv3,tur(i,j,k),i,j,k
                end if
                gv=cb1*q16(i,j,k)*svot
                !********生成项gv考虑旋转和曲率后***********
                if(m==1 .and. n==2)then!初始时不能计算，否则静叶发散
                else
                    s11=gradc(1,i,j,k)
                    s22=gradc(5,i,j,k)
                    s33=gradc(9,i,j,k)
                    s12=0.5d0*(gradc(2,i,j,k)+gradc(4,i,j,k))
                    s13=0.5d0*(gradc(3,i,j,k)+gradc(7,i,j,k))
                    s23=0.5d0*(gradc(6,i,j,k)+gradc(8,i,j,k))
                    ss2=4.d0*(s12*s12+s13*s13+s23*s23)+2.d0*(s11*s11+s22*s22+s33*s33)
                    sss=sqrt(ss2)
                    rp1=rpm(rm)
                    rp2=0.d0
                    rp3=0.d0
                    w12=0.5d0*(gradc(2,i,j,k)-gradc(4,i,j,k))-rp3
                    w13=0.5d0*(gradc(3,i,j,k)-gradc(7,i,j,k))+rp2
                    w23=0.5d0*(gradc(6,i,j,k)-gradc(8,i,j,k))-rp1
                    w21=-w12
                    w31=-w13
                    w32=-w23
                    ww2=4.d0*(w12*w12+w13*w13+w23*w23)
                    ww=sqrt(ww2)               !论文里红圈圈
                    fr1r1=sss/ww  !r星
                    !dd=0.5*(ss2+ww2)  !SA模型文献00年
                    dd=max(ss2,0.09d0*tur(i,j,k)*tur(i,j,k))  !SA模型文献09年
                    fr1r2=2.d0/dd/sqrt(dd)/ww*(gradcs(1,i,j,k)*(w12*s12+w13*s13)+(0.5d0*(gradcs(2,i,j,k)+gradcs(4,i,j,k))-s13*rp1)*(w12*s22+w13*s23+w21*s11+w23*s13)&
                         +(0.5d0*(gradcs(3,i,j,k)+gradcs(7,i,j,k))+s12*rp1)*(w12*s23+w13*s33+w31*s11+w32*s12)+(gradcs(5,i,j,k)-2.d0*s23*rp1)*(w21*s12+w23*s23)&
                         +(0.5d0*(gradcs(6,i,j,k)+gradcs(8,i,j,k))+(s22-s33)*rp1)*(w21*s13+w23*s33+w31*s12+w32*s22)+(gradcs(9,i,j,k)+2.d0*s23*rp1)*(w31*s13+w32*s23))
                    fr1=(1.d0+cr1)*2.d0*fr1r1/(1.d0+fr1r1)*(1.d0-cr3*atan(cr2*fr1r2))-cr1
                    gv=gv*fr1
                end if
            !********耗散项Yv计算***********
                ra=vm/svot
                ga=ra+cw2*(ra**6-ra)
                fw=((ga**(-6)+cw3**(-6))/(1.d0+cw3**(-6)))**(-1.d0/6.d0)
                yv=cw1*fw*q16(i,j,k)*vm*kap*kap
            !********源项计算，只有第6项，所以纳入耗散通量里***********
                qv6(i,j,k)=qv6(i,j,k)+(gv-yv+dv)*vv0(i,j,k)
            end do
        end do
    end do
    deallocate(tur)
    deallocate(pwx)
    deallocate(pwy)
    deallocate(pwz)
    deallocate(gradc)
    deallocate(gradcs)
    end subroutine SAsource
    
subroutine gradsfaceI(si,sj,sk,q,Idqd)
    use global
    implicit none
    real(8),intent(in)  ::si(0:nx+2,1:ny,1:nz),sj(1:nx,0:ny+2,1:nz),sk(1:nx,1:ny,0:nz+2),q(0:nx+1,0:ny+1,0:nz+1)
    real(8),intent(out) ::Idqd(1:nx+1,0:ny+1,0:nz+1)   !和下面J区别是，第一值I由1—nx+1,J由0—nx+1;第二值I从0开始，J从1开始。
    real(8) ::sx,fg,qav,rvol,sx1,sx2,qav1,qav2

    Idqd=0.d0
    do k=1,nz
        do j=1,ny
            do i=1,nx    ! 左右面通量x
                sx= 0.5D0*(si(i,j,k)+si(i+1,j,k))        !此时流过网格中心的通量=当前网格中心i的面积*当前网格的物理量
                fg=q(i,j,k)*sx
            
                Idqd(i,j,k)  = Idqd(i,j,k)  - fg            !此为该网格i,j左界面经过通量以后的值
                Idqd(i+1,j,k)= Idqd(i+1,j,k)+ fg                         !右
            end do
    ! - treat boundary i=1
            sx  =  si(1,j,k)
            fg  =  0.5D0*(q(0,j,k)+q(1,j,k))*sx
            Idqd(1,j,k) = Idqd(1,j,k) + fg
        ! - treat boundary i=il
            sx  =  si(nx+1,j,k)
            fg  =  0.5D0*(q(nx,j,k)+q(nx+1,j,k))*sx
            Idqd(nx+1,j,k) = Idqd(nx+1,j,k) - fg
         end do !i
    end do !j
    
    do k=1,nz              !上下面通量y
        do j=1,ny+1
            do i=2,nx
                sx1=0.5D0*sj(i,j,k)
                qav1=0.5D0*(q(i,j,k)+q(i,j-1,k))
                sx2=0.5D0*sj(i-1,j,k)
                qav2=0.5D0*(q(i-1,j,k)+q(i-1,j-1,k))
                fg=sx1*qav1+sx2*qav2
                Idqd(i,j,k)  = Idqd(i,j,k)  + fg          !?
                Idqd(i,j-1,k)= Idqd(i,j-1,k)- fg
            end do
    ! - treat boundary i=1 最左侧的上下面通量
            sx  = sj(1,j,k)
            qav = 0.5D0*(q(1,j,k)+q(1,j-1,k))
            fg  = qav*sx
            Idqd(1,j,k)   = Idqd(1,j,k)   + fg
            Idqd(1,j-1,k) = Idqd(1,j-1,k) - fg
    ! - treat boundary i=nx+1最右侧的上下面通量
            sx  = sj(nx,j,k)      !此时流过网格中心通量=当前网格中心i与i-1的下侧面积加权*当前网格与其它左、下侧网格的物理量加权
            qav = 0.5D0*(q(nx,j,k)+q(nx,j-1,k))
            fg  = qav*sx
            Idqd(nx+1,j,k)   = Idqd(nx+1,j,k)   + fg
            Idqd(nx+1,j-1,k) = Idqd(nx+1,j-1,k) - fg
        end do
    end do !j
   
    do j=1,ny               !前后面通量z
        do k=1,nz+1
            do i=2,nx
                sx1=0.5D0*sk(i,j,k)
                qav1=0.5D0*(q(i,j,k)+q(i,j,k-1))
                sx2=0.5D0*sk(i-1,j,k)
                qav2=0.5D0*(q(i-1,j,k)+q(i-1,j,k-1))
                fg=sx1*qav1+sx2*qav2
                Idqd(i,j,k)  = Idqd(i,j,k)  + fg          !?
                Idqd(i,j,k-1)= Idqd(i,j,k-1)- fg
            end do
    ! - treat boundary i=1 最左侧的上下面通量
            sx  = sk(1,j,k)
            qav = 0.5D0*(q(1,j,k)+q(1,j,k-1))
            fg  = qav*sx
            Idqd(1,j,k)   = Idqd(1,j,k)   + fg
            Idqd(1,j,k-1) = Idqd(1,j,k-1) - fg
    ! - treat boundary i=nx+1最右侧的上下面通量
            sx  = sk(nx,j,k)      !此时流过网格中心通量=当前网格中心i与i-1的下侧面积加权*当前网格与其它左、下侧网格的物理量加权
            qav = 0.5D0*(q(nx,j,k)+q(nx,j,k-1))
            fg  = qav*sx
            Idqd(nx+1,j,k)   = Idqd(nx+1,j,k)   + fg
            Idqd(nx+1,j,k-1) = Idqd(nx+1,j,k-1) - fg
        end do
    end do !j
    
    do k=1,nz
        do j=1,ny
            do i=2,nx
                rvol       = 2.D0/(vv0(i,j,k)+vv0(i-1,j,k))
                Idqd(i,j,k)= Idqd(i,j,k)*rvol
            end do
            rvol       = 1.D0/vv0(1,j,k)
            Idqd(1,j,k)= Idqd(1,j,k)*rvol
        
            rvol          = 1.D0/vv0(nx,j,k)
            Idqd(nx+1,j,k)= Idqd(nx+1,j,k)*rvol
        end do
    end do
    end subroutine gradsfaceI
    
subroutine gradsfaceJ(si,sj,sk,q,Jdqd)
    use global
    implicit none
    real(8),intent(in)  ::si(0:nx+2,1:ny,1:nz),sj(1:nx,0:ny+2,1:nz),sk(1:nx,1:ny,0:nz+2),q(0:nx+1,0:ny+1,0:nz+1)
    real(8),intent(out) ::Jdqd(0:nx+1,1:ny+1,0:nz+1)   !和下面J区别是，第一值I由1—nx+1,J由0—nx+1;第二值I从0开始，J从1开始。
    real(8) ::sx,fg,qav,rvol,sx1,sx2,qav1,qav2
    
    Jdqd=0.d0
    do k=1,nz
        do i=1,nx
            do j=1,ny
    ! --- bottom face of auxiliary control volume
                sx   =  0.5D0*(sj(i,j,k)+sj(i,j+1,k))
                fg   =  q(i,j,k)*sx
                Jdqd(i,j,k)   = Jdqd(i,j,k)   - fg
                Jdqd(i,j+1,k) = Jdqd(i,j+1,k) + fg
            end do
        ! - treat boundary j=1
            sx  =  sj(i,1,k)
            fg  =  0.5D0*( q(i,0,k)+ q(i,1,k))*sx    !界面值
            Jdqd(i,1,k) =  Jdqd(i,1,k)  +fg
        ! - treat boundary j=jl
            sx  = sj(i,ny+1,k)
            fg  =  0.5D0*( q(i,ny,k)+ q(i,ny+1,k))*sx
            Jdqd(i,ny+1,k) =  Jdqd(i,ny+1,k)  -fg
         end do !i
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
    
    ! divide by the volume of auxiliary cell
    do k=1,nz
        do i=1,nx
            do j=2,ny
                rvol     = 2.D0/(vv0(i,j,k)+vv0(i,j-1,k))
                Jdqd(i,j,k)= Jdqd(i,j,k)*rvol
            end do
            rvol        = 1.D0/vv0(i,1,k)
            Jdqd(i,1,k)   = Jdqd(i,1,k)*rvol
        
            rvol        = 1.D0/vv0(i,ny,k)
            Jdqd(i,ny+1,k)  = Jdqd(i,ny+1,k)*rvol
        end do
    end do
   end subroutine gradsfaceJ
    
subroutine gradsfaceK(si,sj,sk,q,Kdqd)
    use global
    implicit none
    real(8),intent(in)  ::si(0:nx+2,1:ny,1:nz),sj(1:nx,0:ny+2,1:nz),sk(1:nx,1:ny,0:nz+2),q(0:nx+1,0:ny+1,0:nz+1)
    real(8),intent(out) ::Kdqd(0:nx+1,0:ny+1,1:nz+1)   !和下面J区别是，第一值I由1—nx+1,J由0—nx+1;第二值I从0开始，J从1开始。
    real(8) ::sx,fg,qav,rvol,sx1,sx2,qav1,qav2
    
    Kdqd=0.d0
    do j=1,ny
        do i=1,nx
            do k=1,nz
    ! --- right face of auxiliary control volume前后面通量z
                sx= 0.5D0*(sk(i,j,k)+sk(i,j,k+1))        !此时流过网格中心的通量=当前网格中心i的面积*当前网格的物理量
                fg=q(i,j,k)*sx
            
                Kdqd(i,j,k)  = Kdqd(i,j,k)  - fg            !此为该网格i,j左界面经过通量以后的值
                Kdqd(i,j,k+1)= Kdqd(i,j,k+1)+ fg
            end do
    ! - treat boundary i=1
            sx  =  sk(i,j,1)
            fg=0.5*(q(i,j,0)+q(i,j,1))*sx          !固壁上用边界值
            Kdqd(i,j,1) = Kdqd(i,j,1) + fg
        ! - treat boundary i=il
            sx  =  sk(i,j,nz+1)
            fg=0.5*(q(i,j,nz)+q(i,j,nz+1))*sx
            Kdqd(i,j,nz+1) = Kdqd(i,j,nz+1) - fg
         end do !i
    end do !j

    do j=1,ny
        do i=1,nx+1   !左右面通量x
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
                Kdqd(i,j,k)  = Kdqd(i,j,k)  + fg          !?
                Kdqd(i,j-1,k)= Kdqd(i,j-1,k)- fg
            end do
    ! - treat boundary i=1 最左侧的上下面通量
            sx  = sj(i,j,1)
            qav = 0.5D0*(q(i,j,1)+q(i,j-1,1))
            fg  = qav*sx
            Kdqd(i,j,1)   = Kdqd(i,j,1)   + fg
            Kdqd(i,j-1,1) = Kdqd(i,j-1,1) - fg
    ! - treat boundary i=nx+1最右侧的上下面通量
            sx  = sj(i,j,nz)      !此时流过网格中心通量=当前网格中心i与i-1的下侧面积加权*当前网格与其它左、下侧网格的物理量加权
            qav = 0.5D0*(q(i,j,nz)+q(i,j-1,nz))
            fg  = qav*sx
            Kdqd(i,j,nz+1)   = Kdqd(i,j,nz+1)   + fg
            Kdqd(i,j-1,nz+1) = Kdqd(i,j-1,nz+1) - fg
        end do
    end do !j
   
    do j=1,ny
        do i=1,nx
            do k=2,nz
                rvol       = 2.D0/(vv0(i,j,k)+vv0(i,j,k-1))
                Kdqd(i,j,k)= Kdqd(i,j,k)*rvol
            end do
            rvol       = 1.D0/vv0(i,j,1)
            Kdqd(i,j,1)= Kdqd(i,j,1)*rvol
        
            rvol          = 1.D0/vv0(i,j,nz)
            Kdqd(i,j,nz+1)= Kdqd(i,j,nz+1)*rvol
        end do
    end do
    end subroutine gradsfaceK
        
    subroutine gradscentre(direction,q,dqd) !计算网格中心物理量梯度（求湍流粘性项gv的S时用到）
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
                dqd(i,j,k)=dqd(i,j,k)/vv0(i,j,k)
            end do
        end do
    end do
    end subroutine gradscentre
    
subroutine dsdt(q,sss)
    use global
    implicit none
    real(8),intent(in)  ::q(0:nx+1,0:ny+1,0:nz+1)
    real(8),intent(out) ::sss(0:nx+1,0:ny+1,0:nz+1)
    real(8) ::vf,rf,qq1,flu
    
    sss=0.
    !******************x方向********************
    do k=1,nz
        do j=1,ny
            do i=1,nx+1
               vx=0.5d0*(pvx(i,j,k)+pvx(i-1,j,k))
               vy=0.5d0*(pvy(i,j,k)+pvy(i-1,j,k))
               vz=0.5d0*(pvz(i,j,k)+pvz(i-1,j,k))
               vf=vx*s2x(i,j,k)+vy*s2y(i,j,k)+vz*s2z(i,j,k)
               rf=rpm(rm)*zz02(i,j,k)*s2y(i,j,k)-rpm(rm)*yy02(i,j,k)*s2z(i,j,k)     !牵连速度通量
               vf=vf+rf                                                 !相对速度通量
               qq1=0.5d0*(q(i,j,k)+q(i-1,j,k))
               flu=qq1*vf
               sss(i,j,k)=sss(i,j,k)+flu
               sss(i-1,j,k)=sss(i-1,j,k)-flu
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
               rf=rpm(rm)*zz03(i,j,k)*s3y(i,j,k)-rpm(rm)*yy03(i,j,k)*s3z(i,j,k)     !牵连速度通量
               vf=vf+rf                                                 !相对速度通量
               qq1=0.5d0*(q(i,j,k)+q(i,j-1,k))
               flu=qq1*vf
               sss(i,j,k)=sss(i,j,k)+flu
               sss(i,j-1,k)=sss(i,j-1,k)-flu
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
                rf=rpm(rm)*zz01(i,j,k)*s1y(i,j,k)-rpm(rm)*yy01(i,j,k)*s1z(i,j,k)     !牵连速度通量
                vf=vf+rf                                                 !相对速度通量
                qq1=0.5d0*(q(i,j,k)+q(i,j,k-1))
                flu=qq1*vf
                sss(i,j,k)=sss(i,j,k)+flu
                sss(i,j,k-1)=sss(i,j,k-1)-flu
            end do
        end do
    end do
    
    do k=1,nz
        do j=1,ny
            do i=1,nx
                sss(i,j,k)=sss(i,j,k)/vv0(i,j,k)
            end do
        end do
    end do
    end subroutine dsdt
    
subroutine viscosity(temp,q6,cv,kc) !计算粘性系数和导热系数
    use global
    implicit none
    real(8),intent(in)  ::temp,q6
    real(8),intent(out) ::cv,kc
    real(8) ::cvt,fv1,tem
    
    !分子粘性cvl计算
    cvl=cvl0*((temp/t0)**1.5)*(t0+ts)/(temp+ts)
  !湍流粘性cvt计算
    tem=q6/cvl        !此时tem不是温度，而是粘度系数，形状类似于x
    fv1=1.d0/(1.d0+(cv1/tem)**3)
    cvt=q6*fv1
  !粘性系数计算
    cv=cvl+cvt
  !导热系数计算
    kc=cp*(cvl/prl+cvt/prt)
    end subroutine viscosity
    
subroutine pred(ims) !每步R-K推进方法
    use global
    implicit none
    integer,intent(in) ::ims  !ims=1,奇数步进行隐式残差光顺
    real(8),allocatable ::py1(:,:,:),py2(:,:,:),py3(:,:,:),py4(:,:,:),py5(:,:,:),py6(:,:,:)
    real(8),allocatable ::ax(:),bx(:),ay(:),by(:),az(:),bz(:),cx(:),dx(:),cy(:),dy(:),cz(:),dz(:)

    allocate(py1(nx,ny,nz))
    allocate(py2(nx,ny,nz))
    allocate(py3(nx,ny,nz))
    allocate(py4(nx,ny,nz))
    allocate(py5(nx,ny,nz))
    allocate(py6(nx,ny,nz))
    do k=1,nz!RK待加项即delU = -a*deltao*R'                     !计算方程右边项括号内容R'=a3*u+ts+delt*R）
        do j=1,ny
            do i=1,nx
                py1(i,j,k)=-timl*time(i,j,k)/(delt+aa3*time(i,j,k))*(aa3*q01(i,j,k)+ts1(i,j,k)+(qc1(i,j,k)-av1(i,j,k))*delt/vv0(i,j,k))
                py2(i,j,k)=-timl*time(i,j,k)/(delt+aa3*time(i,j,k))*(aa3*q02(i,j,k)+ts2(i,j,k)+(qc2(i,j,k)-av2(i,j,k)-qv2(i,j,k))*delt/vv0(i,j,k))
                py3(i,j,k)=-timl*time(i,j,k)/(delt+aa3*time(i,j,k))*(aa3*q03(i,j,k)+ts3(i,j,k)+(qc3(i,j,k)-av3(i,j,k)-qv3(i,j,k))*delt/vv0(i,j,k))
                py4(i,j,k)=-timl*time(i,j,k)/(delt+aa3*time(i,j,k))*(aa3*q04(i,j,k)+ts4(i,j,k)+(qc4(i,j,k)-av4(i,j,k)-qv4(i,j,k))*delt/vv0(i,j,k))
                py5(i,j,k)=-timl*time(i,j,k)/(delt+aa3*time(i,j,k))*(aa3*q05(i,j,k)+ts5(i,j,k)+(qc5(i,j,k)-av5(i,j,k)-qv5(i,j,k))*delt/vv0(i,j,k))
                py6(i,j,k)=-timl*time(i,j,k)/(delt+aa3*time(i,j,k))*(aa3*q06(i,j,k)+ts6(i,j,k)+(qc6(i,j,k)-av6(i,j,k)-qv6(i,j,k))*delt/vv0(i,j,k))
            end do
        end do
    end do
    ta=1.d0
    if(ims==1)then!隐式残差光顺ave
        allocate(ax(nx))
        allocate(bx(nx))
        allocate(ay(ny))
        allocate(by(ny))
        allocate(az(nz))
        allocate(bz(nz))
        allocate(cx(nx))
        allocate(dx(nx))
        allocate(cy(ny))
        allocate(dy(ny))
        allocate(cz(nz))
        allocate(dz(nz))
        ax=-ta
        ay=-ta
        az=-ta
        bx=1.d0+2.d0*ta
        by=1.d0+2.d0*ta
        bz=1.d0+2.d0*ta
        cx=-1.5*ta
        cy=-1.5*ta
        cz=-1.5*ta
        dx=1.d0+3.d0*ta
        dy=1.d0+3.d0*ta
        dz=1.d0+3.d0*ta
        do k=1,nz                                            !x方向光顺
            do j=1,ny
                call tdma(ax(2:nx),bx(1:nx),ax(1:nx-1),py1(1:nx,j,k),py1(1:nx,j,k),nx)
                call tdma(ax(2:nx),bx(1:nx),ax(1:nx-1),py2(1:nx,j,k),py2(1:nx,j,k),nx)
                call tdma(ax(2:nx),bx(1:nx),ax(1:nx-1),py3(1:nx,j,k),py3(1:nx,j,k),nx)
                call tdma(ax(2:nx),bx(1:nx),ax(1:nx-1),py4(1:nx,j,k),py4(1:nx,j,k),nx)
                call tdma(ax(2:nx),bx(1:nx),ax(1:nx-1),py5(1:nx,j,k),py5(1:nx,j,k),nx)
                call tdma(cx(2:nx),dx(1:nx),cx(1:nx-1),py6(1:nx,j,k),py6(1:nx,j,k),nx)
            end do
        end do
        do k=1,nz                                            !y方向光顺
            do i=1,nx
                call tdma(ay(2:ny),by(1:ny),ay(1:ny-1),py1(i,1:ny,k),py1(i,1:ny,k),ny)
                call tdma(ay(2:ny),by(1:ny),ay(1:ny-1),py2(i,1:ny,k),py2(i,1:ny,k),ny)
                call tdma(ay(2:ny),by(1:ny),ay(1:ny-1),py3(i,1:ny,k),py3(i,1:ny,k),ny)
                call tdma(ay(2:ny),by(1:ny),ay(1:ny-1),py4(i,1:ny,k),py4(i,1:ny,k),ny)
                call tdma(ay(2:ny),by(1:ny),ay(1:ny-1),py5(i,1:ny,k),py5(i,1:ny,k),ny)
                call tdma(cy(2:ny),dy(1:ny),cy(1:ny-1),py6(i,1:ny,k),py6(i,1:ny,k),ny)
            end do
        end do
        do j=1,ny                                            !z方向光顺
            do i=1,nx
                call tdma(az(2:nz),bz(1:nz),az(1:nz-1),py1(i,j,1:nz),py1(i,j,1:nz),nz)
                call tdma(az(2:nz),bz(1:nz),az(1:nz-1),py2(i,j,1:nz),py2(i,j,1:nz),nz)
                call tdma(az(2:nz),bz(1:nz),az(1:nz-1),py3(i,j,1:nz),py3(i,j,1:nz),nz)
                call tdma(az(2:nz),bz(1:nz),az(1:nz-1),py4(i,j,1:nz),py4(i,j,1:nz),nz)
                call tdma(az(2:nz),bz(1:nz),az(1:nz-1),py5(i,j,1:nz),py5(i,j,1:nz),nz)
                call tdma(cz(2:nz),dz(1:nz),cz(1:nz-1),py6(i,j,1:nz),py6(i,j,1:nz),nz)
            end do
        end do
        deallocate(ax)
        deallocate(bx)
        deallocate(ay)
        deallocate(by)
        deallocate(az)
        deallocate(bz)
        deallocate(cx)
        deallocate(dx)
        deallocate(cy)
        deallocate(dy)
        deallocate(cz)
        deallocate(dz)
    end if
    q11(1:nx,1:ny,1:nz)=q01(1:nx,1:ny,1:nz)+py1(1:nx,1:ny,1:nz)
    q12(1:nx,1:ny,1:nz)=q02(1:nx,1:ny,1:nz)+py2(1:nx,1:ny,1:nz)
    q13(1:nx,1:ny,1:nz)=q03(1:nx,1:ny,1:nz)+py3(1:nx,1:ny,1:nz)
    q14(1:nx,1:ny,1:nz)=q04(1:nx,1:ny,1:nz)+py4(1:nx,1:ny,1:nz)
    q15(1:nx,1:ny,1:nz)=q05(1:nx,1:ny,1:nz)+py5(1:nx,1:ny,1:nz)
    q16(1:nx,1:ny,1:nz)=q06(1:nx,1:ny,1:nz)+py6(1:nx,1:ny,1:nz)
    deallocate(py1)
    deallocate(py2)
    deallocate(py3)
    deallocate(py4)
    deallocate(py5)
    deallocate(py6)
    end subroutine pred
    
subroutine tdma(a,b,c,d,x,n) !TDMA算法
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
    
    subroutine residual
    use global
    implicit none
    
    rmsm=0.d0
    do k=1,nz
        do j=1,ny
            do i=1,nx
                rmsm=rmsm+((Q11(i,j,k)-Q01(i,j,k))/time(i,j,k))**2
            end do
        end do
    end do
    rmsm=0.5d0*log10(rmsm/dble(nx*ny*nz))
end subroutine residual

subroutine store!存储前两步真实时步值，用于真实时间导数项
    use global
    implicit none

    q21(1:nx,1:ny,1:nz,1)=q21(1:nx,1:ny,1:nz,2)
    q22(1:nx,1:ny,1:nz,1)=q22(1:nx,1:ny,1:nz,2)
    q23(1:nx,1:ny,1:nz,1)=q23(1:nx,1:ny,1:nz,2)
    q24(1:nx,1:ny,1:nz,1)=q24(1:nx,1:ny,1:nz,2)
    q25(1:nx,1:ny,1:nz,1)=q25(1:nx,1:ny,1:nz,2)
    q26(1:nx,1:ny,1:nz,1)=q26(1:nx,1:ny,1:nz,2)
    
    q21(1:nx,1:ny,1:nz,2)=q11(1:nx,1:ny,1:nz)
    q22(1:nx,1:ny,1:nz,2)=q12(1:nx,1:ny,1:nz)
    q23(1:nx,1:ny,1:nz,2)=q13(1:nx,1:ny,1:nz)
    q24(1:nx,1:ny,1:nz,2)=q14(1:nx,1:ny,1:nz)
    q25(1:nx,1:ny,1:nz,2)=q15(1:nx,1:ny,1:nz)
    q26(1:nx,1:ny,1:nz,2)=q16(1:nx,1:ny,1:nz)
    
    open(myid+300,file='pause-'//trim(adjustl(id_m))//'myid.dat',form="unformatted")!记录下文件，方便中断后再读入
    rewind(myid+300) !覆盖之前的文件
    write(myid+300) m
    write(myid+300) n
    write(myid+300) nitt
    write(myid+300) q21(1:nx,1:ny,1:nz,1:2)
    write(myid+300) q22(1:nx,1:ny,1:nz,1:2)
    write(myid+300) q23(1:nx,1:ny,1:nz,1:2)
    write(myid+300) q24(1:nx,1:ny,1:nz,1:2)
    write(myid+300) q25(1:nx,1:ny,1:nz,1:2)
    write(myid+300) q26(1:nx,1:ny,1:nz,1:2)
    close(myid+300)
    end subroutine store

subroutine probe
    use global
    implicit none
    
    j=ny-2
    k=5
    if(rm==1)then!动叶，静叶前后缘每个叶道顶部都布置探针
        if(myid==3+lbm*numpp)then
            i=nx-1
            vx=q12(i,j,k)/q11(i,j,k)
            vy=q13(i,j,k)/q11(i,j,k)
            vz=q14(i,j,k)/q11(i,j,k)
            qq2=vx*vx+vy*vy+vz*vz
            pp=0.4d0*(q15(i,j,k)-0.5d0*q11(i,j,k)*qq2)   !静压
            write(131,"(i10,F15.6)")nitt,vx
            write(132,"(i10,F15.6)")nitt,pp
        end if
        if(myid==9+lbm*numpp)then
            i=2
            vx=q12(i,j,k)/q11(i,j,k)
            vy=q13(i,j,k)/q11(i,j,k)
            vz=q14(i,j,k)/q11(i,j,k)
            qq2=vx*vx+vy*vy+vz*vz
            pp=0.4d0*(q15(i,j,k)-0.5d0*q11(i,j,k)*qq2)   !静压
            write(133,"(i10,F15.6)")nitt,vx
            write(134,"(i10,F15.6)")nitt,pp
        end if
    end if
    if(rm==2)then
        if(myid==1+lbm*numpp)then
            i=nx-1
            vx=q12(i,j,k)/q11(i,j,k)
            vy=q13(i,j,k)/q11(i,j,k)
            vz=q14(i,j,k)/q11(i,j,k)
            qq2=vx*vx+vy*vy+vz*vz
            pp=0.4d0*(q15(i,j,k)-0.5d0*q11(i,j,k)*qq2)   !静压
            write(231,"(i10,F15.6)")nitt,vx
            write(232,"(i10,F15.6)")nitt,pp
        end if
        if(myid==8+lbm*numpp)then
            i=2
            vx=q12(i,j,k)/q11(i,j,k)
            vy=q13(i,j,k)/q11(i,j,k)
            vz=q14(i,j,k)/q11(i,j,k)
            qq2=vx*vx+vy*vy+vz*vz
            pp=0.4d0*(q15(i,j,k)-0.5d0*q11(i,j,k)*qq2)   !静压
            write(233,"(i10,F15.6)")nitt,vx
            write(234,"(i10,F15.6)")nitt,pp
        end if
    end if
    end subroutine probe

subroutine test
    use global
    implicit none
    real(8) :: temp,t1,t2,t3,t4,t5,t6,ptur,vf,sx,sr,sxx,srr,s1,s2,s3,vr,va,dx,dy,dz,dr,tn,tt,vn,vt,ss,ht1,ht0,flu1,flu2,flu3,flu4,flu5,flu6,flu11,flu22,flu33,flu44,flu55,flu66,flu01,flu02,flu03,flu04,flu05,flu06
    real(8),allocatable :: test01(:),test02(:),test03(:),test04(:),test05(:),test06(:)
    integer :: myidll
    
    if(slidm==1)then
        i=nx
        ii=nx+1
    else
        i=0
        ii=1
    end if
    flu1=0.
    flu2=0.
    flu3=0.
    flu4=0.
    flu5=0.
    flu6=0.
    do j=1,ny
        temp=0.
        t1=0.
        t2=0.
        t3=0.
        t4=0.
        t5=0.
        t6=0.
        do k=1,nz
            y1=yy0(i,j,k)
            z1=zz0(i,j,k)
            rr=sqrt(y1*y1+z1*z1)
            sir=z1/rr
            cor=y1/rr
            dim=Q11(i,j,k)
            vx=q12(i,j,k)/Q11(i,j,k)
            vy=q13(i,j,k)/Q11(i,j,k)
            vz=q14(i,j,k)/Q11(i,j,k)
            en=Q15(i,j,k)
            ptur=Q16(i,j,k)
            va=vz*cor-vy*sir
            vr=vz*sir+vy*cor
            qq2=vx*vx+vy*vy+vz*vz
            pp=0.4*(en-0.5*dim*qq2)
            ht0=en+pp
            s1=-0.5*(s2x(ii,j,k)+s2x(i,j,k))
            s2=-0.5*(s2y(ii,j,k)+s2y(i,j,k))
            s3=-0.5*(s2z(ii,j,k)+s2z(i,j,k))
            ss=sqrt(s1*s1+s2*s2+s3*s3)
            dx=s1/ss
            dy=s2/ss
            dz=s3/ss
            dr=sqrt(1-dx*dx)
            vn=vx*dx+vr*dr
            vt=-vx*dr+vr*dx
            vf=vx*s1+vy*s2+vz*s3
            
            temp=temp+ss
            t1=t1+dim*vn*ss    !类似于对流通量的各项
            t2=t2+(dim*vn*vx+pp*dx)*ss
            t3=t3+dim*vn*va*ss  !周向
            t4=t4+(dim*vn*vr+pp*dr)*ss
            t5=t5+vn*ht0*ss
            t6=t6+vn*ptur*ss
        end do
        flu1=flu1+t1/temp
        flu2=flu2+t2/temp
        flu3=flu3+t3/temp
        flu4=flu4+t4/temp
        flu5=flu5+t5/temp
        flu6=flu6+t6/temp
    end do
    !处理总通量
    allocate(test01(0:lb(1)+lb(2)-1))
    allocate(test02(0:lb(1)+lb(2)-1))
    allocate(test03(0:lb(1)+lb(2)-1))
    allocate(test04(0:lb(1)+lb(2)-1))
    allocate(test05(0:lb(1)+lb(2)-1))
    allocate(test06(0:lb(1)+lb(2)-1))
    if(rm==1)then
        myidll=yl+xln(0)*yln+xln(1)*yln
    else if(rm==2)then
        myidll=yl+lb(rm-1)*numpp
    end if
    myidl=myidll+numpp
    myidr=myidll+(lb(rm)-1)*numpp
    if(myid==myidll)then
        j=myid/numpp
        test01(j)=flu1
        test02(j)=flu2
        test03(j)=flu3
        test04(j)=flu4
        test05(j)=flu5
        test06(j)=flu6
        do i=myidl,myidr,numpp
            j=i/numpp
            call MPI_RECV(test01(j),1,MPI_DOUBLE_PRECISION,i,43,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(test02(j),1,MPI_DOUBLE_PRECISION,i,44,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(test03(j),1,MPI_DOUBLE_PRECISION,i,45,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(test04(j),1,MPI_DOUBLE_PRECISION,i,46,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(test05(j),1,MPI_DOUBLE_PRECISION,i,47,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(test06(j),1,MPI_DOUBLE_PRECISION,i,48,MPI_COMM_WORLD,status,ierr)
        end do
        flu11=0.
        flu22=0.
        flu33=0.
        flu44=0.
        flu55=0.
        flu66=0.
        do i=myidll,myidr,numpp
            j=i/numpp
            flu11=flu11+test01(j)
            flu22=flu22+test02(j)
            flu33=flu33+test03(j)
            flu44=flu44+test04(j)
            flu55=flu55+test05(j)
            flu66=flu66+test06(j)
        end do
    else
        call MPI_SEND(flu1,1,MPI_DOUBLE_PRECISION,myidll,43,MPI_COMM_WORLD,ierr)
        call MPI_SEND(flu2,1,MPI_DOUBLE_PRECISION,myidll,44,MPI_COMM_WORLD,ierr)
        call MPI_SEND(flu3,1,MPI_DOUBLE_PRECISION,myidll,45,MPI_COMM_WORLD,ierr)
        call MPI_SEND(flu4,1,MPI_DOUBLE_PRECISION,myidll,46,MPI_COMM_WORLD,ierr)
        call MPI_SEND(flu5,1,MPI_DOUBLE_PRECISION,myidll,47,MPI_COMM_WORLD,ierr)
        call MPI_SEND(flu6,1,MPI_DOUBLE_PRECISION,myidll,48,MPI_COMM_WORLD,ierr)
    end if
    if(lm==0)then
        if(yl==0)then
            call MPI_RECV(flu01,1,MPI_DOUBLE_PRECISION,myid+xln(xlln),53,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(flu02,1,MPI_DOUBLE_PRECISION,myid+xln(xlln),54,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(flu03,1,MPI_DOUBLE_PRECISION,myid+xln(xlln),55,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(flu04,1,MPI_DOUBLE_PRECISION,myid+xln(xlln),56,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(flu05,1,MPI_DOUBLE_PRECISION,myid+xln(xlln),57,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(flu06,1,MPI_DOUBLE_PRECISION,myid+xln(xlln),58,MPI_COMM_WORLD,status,ierr)
            write(121,*)nitt,(flu11+flu01)/lb(rm)
            write(122,*)nitt,(flu22+flu02)/lb(rm)
            write(123,*)nitt,(flu33+flu03)/lb(rm)
            write(124,*)nitt,(flu44+flu04)/lb(rm)
            write(125,*)nitt,(flu55+flu05)/lb(rm)
            write(126,*)nitt,(flu66+flu06)/lb(rm)
        else
            call MPI_SEND(flu11,1,MPI_DOUBLE_PRECISION,myid-xln(xlln),53,MPI_COMM_WORLD,ierr)
            call MPI_SEND(flu22,1,MPI_DOUBLE_PRECISION,myid-xln(xlln),54,MPI_COMM_WORLD,ierr)
            call MPI_SEND(flu33,1,MPI_DOUBLE_PRECISION,myid-xln(xlln),55,MPI_COMM_WORLD,ierr)
            call MPI_SEND(flu44,1,MPI_DOUBLE_PRECISION,myid-xln(xlln),56,MPI_COMM_WORLD,ierr)
            call MPI_SEND(flu55,1,MPI_DOUBLE_PRECISION,myid-xln(xlln),57,MPI_COMM_WORLD,ierr)
            call MPI_SEND(flu66,1,MPI_DOUBLE_PRECISION,myid-xln(xlln),58,MPI_COMM_WORLD,ierr)
        end if
    end if
    deallocate(test01)
    deallocate(test02)
    deallocate(test03)
    deallocate(test04)
    deallocate(test05)
    deallocate(test06)
    end subroutine test

 subroutine flow  !进出口质量流量
    use global
    implicit none
    real(8) :: vf1,vf2,qin,qout,qinn,qoutt
    real(8),allocatable :: qin0(:),qout0(:)
    
    if(rm==1 .and. xll==0)then
        qin=0.d0
        do j=1,ny
            do k=1,nz
                vf1=-(s2x(1,j,k)*q12(1,j,k)+s2y(1,j,k)*q13(1,j,k)+s2z(1,j,k)*q14(1,j,k))
                qin=qin+vf1
            end do
        end do
        if(myid==0)then
            allocate(qin0(0:lb(rm)-1))
            k=myid/numpp
            qin0(k)=qin
            qinn=0.0
            do j=myid+numpp,myid+(lb(rm)-1)*numpp,numpp
                k=j/numpp
                call MPI_RECV(qin0(k),1,MPI_DOUBLE_PRECISION,j,9,MPI_COMM_WORLD,status,ierr)
            end do
            do k=0,lb(rm)-1
                qinn=qinn+qin0(k)
            end do
            do j=myid+xln(xlln),myid+xln(xlln)+(lb(rm)-1)*numpp,numpp
                k=j/numpp
                call MPI_RECV(qin0(k),1,MPI_DOUBLE_PRECISION,j,9,MPI_COMM_WORLD,status,ierr)
            end do
            do k=0,lb(rm)-1
                qinn=qinn+qin0(k)
            end do
            deallocate(qin0)
            write(4,"(i5,F15.6)")m,qinn      !进出口质量流量随叶片通过周期的变化
        else
            call MPI_SEND(qin,1,MPI_DOUBLE_PRECISION,0,9,MPI_COMM_WORLD,ierr)
        end if
    end if
    if(rm==2 .and. xll==xln(0)+xln(1)+xln(2)-1)then
        qout=0.d0
        do j=1,ny
            do k=1,nz
                vf2=-(s2x(nx+1,j,k)*q12(nx,j,k)+s2y(nx+1,j,k)*q13(nx,j,k)+s2z(nx+1,j,k)*q14(nx,j,k))
                qout=qout+vf2
            end do
        end do
        if(myid==lb(1)*numpp+xln(0)*yln+xln(1)*yln+xln(2)-1)then
            allocate(qout0(lb(rm-1):lb(rm-1)+lb(rm)-1))
            k=myid/numpp
            qout0(k)=qout
            qoutt=0.0
            do j=myid+numpp,myid+(lb(rm)-1)*numpp,numpp
                k=j/numpp
                call MPI_RECV(qout0(k),1,MPI_DOUBLE_PRECISION,j,19,MPI_COMM_WORLD,status,ierr)
            end do
            do k=lb(rm-1),lb(rm-1)+lb(rm)-1
                qoutt=qoutt+qout0(k)
            end do
            do j=myid+xln(xlln),myid+xln(xlln)+(lb(rm)-1)*numpp,numpp
                k=j/numpp
                call MPI_RECV(qout0(k),1,MPI_DOUBLE_PRECISION,j,19,MPI_COMM_WORLD,status,ierr)
            end do
            do k=lb(rm-1),lb(rm-1)+lb(rm)-1
                qoutt=qoutt+qout0(k)
            end do
            deallocate(qout0)
            write(7,"(i5,F15.6)")m,qoutt      !进出口质量流量随叶片通过周期的变化
        else
            call MPI_SEND(qout,1,MPI_DOUBLE_PRECISION,lb(1)*numpp+xln(0)*yln+xln(1)*yln+xln(2)-1,19,MPI_COMM_WORLD,ierr)
        end if
    end if
    end subroutine flow
    
subroutine in1out !进出口压比，温比，效率
    use global
    implicit none
    integer :: mml
    real(8) :: vf2,h00,t2,p2,d1,d2,d0,h0,h1,h2,pout,tout,d11,d22,eff,sd110,seff00
    real(8),allocatable ::d10(:),eff0(:),d110(:),eff00(:),dd1(:),eeff(:)
    
    if(rm==2 .and. xll==xln(0)+xln(1)+xln(2)-1)then
        allocate(d110(1:ny))
        allocate(eff00(1:ny))
        h00=0.
        h1=0.
        h2=0.
        do j=1,ny
            d0=0.
            d1=0.
            d2=0.
            do k=1,nz
                vf2=-(s2x(nx+1,j,k)*q12(nx,j,k)+s2y(nx+1,j,k)*q13(nx,j,k)+s2z(nx+1,j,k)*q14(nx,j,k))
                vx=q12(nx,j,k)/q11(nx,j,k)
                vy=q13(nx,j,k)/q11(nx,j,k)
                vz=q14(nx,j,k)/q11(nx,j,k)
                qq2=vx*vx+vy*vy+vz*vz
                pp=0.4d0*(q15(nx,j,k)-0.5d0*q11(nx,j,k)*qq2)   !出口静压
                h0=(q15(nx,j,k)+pp)/q11(nx,j,k)    !出口总焓
                t2=h0/cp                           !出口总温
                p2=pp*(1.d0-0.5d0*qq2/h0)**(-3.5)  !出口总压
                d0=d0+vf2
                d1=d1+t2*vf2
                d2=d2+p2*vf2
                h00=h00+vf2
                h1=h1+t2*vf2
                h2=h2+p2*vf2
            end do
            !******气动参数沿叶高分布
            tout=d1/d0
            pout=d2/d0
            d11=pout/pt      !总压比
            d22=tout/(ht/cp) !总温比
            eff=(d11**(2.d0/7.d0)-1.d0)/(d22-1.d0)*100.d0!效率
            if(myid==lb(1)*numpp+xln(0)*yln+xln(1)*yln+(yl+1)*xln(2)-1)then
                allocate(d10(lb(rm-1):lb(rm-1)+lb(rm)-1))
                allocate(eff0(lb(rm-1):lb(rm-1)+lb(rm)-1))
                k=myid/numpp
                d10(k)=d11
                eff0(k)=eff
                do jj=myid+numpp,myid+(lb(rm)-1)*numpp,numpp
                    k=jj/numpp
                    call MPI_RECV(d10(k),1,MPI_DOUBLE_PRECISION,jj,29,MPI_COMM_WORLD,status,ierr)
                    call MPI_RECV(eff0(k),1,MPI_DOUBLE_PRECISION,jj,39,MPI_COMM_WORLD,status,ierr)
                end do
                do k=lb(rm-1),lb(rm-1)+lb(rm)-1
                    d110(j)=d110(j)+d10(k)
                    eff00(j)=eff00(j)+eff0(k)
                end do
                deallocate(d10)
                deallocate(eff0)
                d110(j)=d110(j)/dble(lb(rm))
                eff00(j)=eff00(j)/dble(lb(rm))
            else
                call MPI_SEND(d11,1,MPI_DOUBLE_PRECISION,lb(1)*numpp+xln(0)*yln+xln(1)*yln+(yl+1)*xln(2)-1,29,MPI_COMM_WORLD,ierr)
                call MPI_SEND(eff,1,MPI_DOUBLE_PRECISION,lb(1)*numpp+xln(0)*yln+xln(1)*yln+(yl+1)*xln(2)-1,39,MPI_COMM_WORLD,ierr)
            end if
        end do
        if(lm==0)then
            if(yl==0)then
                allocate(dd1(1:ny0))
                allocate(eeff(1:ny0))
                dd1(1:ny)=d110(1:ny)
                eeff(1:ny)=eff00(1:ny)
                call MPI_RECV(dd1(ny+1:ny0),ny,MPI_DOUBLE_PRECISION,myid+2,49,MPI_COMM_WORLD,status,ierr)
                call MPI_RECV(eeff(ny+1:ny0),ny,MPI_DOUBLE_PRECISION,myid+2,59,MPI_COMM_WORLD,status,ierr)
                open(88,file=trim(adjustl(id_mm))//'period-eff.dat')
                open(89,file=trim(adjustl(id_mm))//'period-pbi.dat')
                do j=1,ny0
                    write(88,"(2F15.6)")eeff(j)/100.,temp0(j)
                    write(89,"(2F15.6)")dd1(j),temp0(j)
                end do
                close(88)
                close(89)
                deallocate(dd1)
                deallocate(eeff)
            else
                call MPI_SEND(d110(1:ny),ny,MPI_DOUBLE_PRECISION,myid-2,49,MPI_COMM_WORLD,ierr)
                call MPI_SEND(eff00(1:ny),ny,MPI_DOUBLE_PRECISION,myid-2,59,MPI_COMM_WORLD,ierr)
            end if
        end if
        deallocate(d110)
        deallocate(eff00)
        !*******整体性能曲线
        tout=h1/h00
        pout=h2/h00
        d11=pout/pt      !总压比
        d22=tout/(ht/cp) !总温比
        eff=(d11**(2.d0/7.d0)-1.d0)/(d22-1.d0)*100.d0!效率
        if(myid==lb(1)*numpp+xln(0)*yln+xln(1)*yln+xln(2)-1)then
            allocate(d10(lb(rm-1):lb(rm-1)+lb(rm)-1))
            allocate(eff0(lb(rm-1):lb(rm-1)+lb(rm)-1))
            k=myid/numpp
            d10(k)=d11
            eff0(k)=eff
            sd110=0.
            seff00=0.
            do j=myid+numpp,myid+(lb(rm)-1)*numpp,numpp
                k=j/numpp
                call MPI_RECV(d10(k),1,MPI_DOUBLE_PRECISION,j,69,MPI_COMM_WORLD,status,ierr)
                call MPI_RECV(eff0(k),1,MPI_DOUBLE_PRECISION,j,79,MPI_COMM_WORLD,status,ierr)
            end do
            do k=lb(rm-1),lb(rm-1)+lb(rm)-1
                sd110=sd110+d10(k)
                seff00=seff00+eff0(k)
            end do
            do j=myid+xln(xlln),myid+xln(xlln)+(lb(rm)-1)*numpp,numpp
                k=j/numpp
                call MPI_RECV(d10(k),1,MPI_DOUBLE_PRECISION,j,69,MPI_COMM_WORLD,status,ierr)
                call MPI_RECV(eff0(k),1,MPI_DOUBLE_PRECISION,j,79,MPI_COMM_WORLD,status,ierr)
            end do
            do k=lb(rm-1),lb(rm-1)+lb(rm)-1
                sd110=sd110+d10(k)
                seff00=seff00+eff0(k)
            end do
            deallocate(d10)
            deallocate(eff0)
            write(8,"(i10,F15.6)")m,sd110/(lb(rm)*yln)
            write(11,"(i10,F15.6)")m,seff00/(lb(rm)*yln)/100.   !进出口面上的物理量随叶片通过周期的变化
        else
            call MPI_SEND(d11,1,MPI_DOUBLE_PRECISION,lb(1)*numpp+xln(0)*yln+xln(1)*yln+xln(2)-1,69,MPI_COMM_WORLD,ierr)
            call MPI_SEND(eff,1,MPI_DOUBLE_PRECISION,lb(1)*numpp+xln(0)*yln+xln(1)*yln+xln(2)-1,79,MPI_COMM_WORLD,ierr)
        end if
    end if
    end subroutine in1out
    
subroutine output0  !输出第nng层的质量流量和各个时间点的流场信息。
    use global
    implicit none
    integer :: chd
       
    write(id_mm,'(i5)') m
    call flow
    call in1out        !求出口参数
    call span(spa1,spa2)         !spa1和spa2叶高截面物理量
    if(rm==1 .and. xll==1)then
        i=nx
        write(secname,"(a2)") 'rl'
        call chord
    end if
    if(rm==1 .and. xll==4)then
        i=2
        write(secname,"(a2)") 'rt'
        call chord
    end if
    if(rm==2 .and. xll==0)then
        i=nx
        write(secname,"(a2)") 'sl'
        call chord
    end if
    if(rm==2 .and. xll==3)then
        i=2
        write(secname,"(a2)") 'st'
        call chord
    end if
    end subroutine output0

subroutine span(spa11,spa22) !该书写只适用于分别大于小于50情况
    use global
    implicit none
    real(8),intent(in) ::spa11,spa22
    real(8) :: temp,qqw,qqr,spa,tem
    real(8),allocatable :: hx(:),hr(:),sss(:),hs(:),hu(:),hp(:),hwma(:),hmiu(:),huu(:,:),hpp(:,:),hwwma(:,:),hss(:,:),hmmiu(:,:)
    
    if(spa11<50.)then!只是个粗略判断标准
        l=0
    else
        l=1
    end if
    if(spa22<50.)then!只是个粗略判断标准
        kk=0
    else
        kk=1
    end if
    allocate(sss(1:ny))
    allocate(hx(1:ny+1))
    allocate(hr(1:ny+1))
    allocate(hu(1:ny))
    allocate(hp(1:ny))
    allocate(hwma(1:ny))
    allocate(hs(1:ny))
    allocate(hmiu(1:ny))
    allocate(huu(1:nx,1:nz))
    allocate(hpp(1:nx,1:nz))
    allocate(hwwma(1:nx,1:nz))
    allocate(hss(1:nx,1:nz))
    allocate(hmmiu(1:nx,1:nz))
    if(yl==l)then
        hxx(1:nx+1,1:nz+1)=hxx1(1:nx+1,1:nz+1)
        hyy(1:nx+1,1:nz+1)=hyy1(1:nx+1,1:nz+1)
        hzz(1:nx+1,1:nz+1)=hzz1(1:nx+1,1:nz+1)
        spa0(1:nx,1:nz)=spa01(1:nx,1:nz)
        spa=spa11
    end if
    if(yl==kk)then
        hxx(1:nx+1,1:nz+1)=hxx2(1:nx+1,1:nz+1)
        hyy(1:nx+1,1:nz+1)=hyy2(1:nx+1,1:nz+1)
        hzz(1:nx+1,1:nz+1)=hzz2(1:nx+1,1:nz+1)
        spa0(1:nx,1:nz)=spa02(1:nx,1:nz)
        spa=spa22
    end if
        
    do k=1,nz
        do i=1,nx
            sss(1)=yl*hspa(i,k)!yl=0进程此处为0，yl=1时起点为原来的ny+1值
            do j=1,ny
                hx(j)=xx0(i,j,k)
                y1=yy0(i,j,k)
                z1=zz0(i,j,k)
                hr(j)=sqrt(y1*y1+z1*z1)
                if(j>1)then
                    sss(j)=sss(j-1)+sqrt((hx(j)-hx(j-1))**2+(hr(j)-hr(j-1))**2) !以上为计算坐标各网格中心坐标
                end if

                hu(j)=Q12(i,j,k)/Q11(i,j,k)
                vy=Q13(i,j,k)/Q11(i,j,k)
                vz=Q14(i,j,k)/Q11(i,j,k)
                qq2=hu(j)*hu(j)+vy*vy+vz*vz!绝对速度平方
                wx=hu(j)
                wy=vy+rpm(rm)*z1
                wz=vz-rpm(rm)*y1
                qqw=wx*wx+wy*wy+wz*wz          !相对速度平方
                hp(j)=0.4d0*(Q15(i,j,k)-0.5d0*Q11(i,j,k)*qq2)!静压
                a=1.4d0*hp(j)/Q11(i,j,k)
                hwma(j)=sqrt(qqw/a)  !相对马赫数
                hs(j)=ccv*log((hp(j)/ppin)/((Q11(i,j,k)/dimin)**1.4))
                tem=hp(j)/(Q11(i,j,k)*rg)!静温
                cvl=cvl0*((tem/t0)**1.5)*(t0+ts)/(tem+ts)!分子粘度
                hmiu(j)=Q16(i,j,k)/cvl           !涡粘度/分子粘度
            end do
            call wl(sss(1:ny),hu(1:ny),spa0(i,k),huu(i,k),ny,1)
            call wl(sss(1:ny),hwma(1:ny),spa0(i,k),hwwma(i,k),ny,1)
            call wl(sss(1:ny),hp(1:ny),spa0(i,k),hpp(i,k),ny,1)
            call wl(sss(1:ny),hs(1:ny),spa0(i,k),hss(i,k),ny,1)
            call wl(sss(1:ny),hmiu(1:ny),spa0(i,k),hmmiu(i,k),ny,1)
        end do
    end do
    !********周期数-90%span-进程号 叶排.dat
    write(nnspan,'(f5.1)') spa
    open(31,file=trim(adjustl(id_mm))//'period-'//trim(adjustl(nnspan))//'%span-'//trim(adjustl(id_m))//'myid.dat')
    write(31,*) "VARIABLES=x,y,z,pvx,wma,pressure,entropy,viscosity"  !抬头
    write(31,"(1x,A,i4,2x,A,i4,2x,A,i4,2x,A)") "ZONE I=",nx+1,"J=",1,"K=",nz+1,"DATAPACKING=BLOCK, VARLOCATION=([4-8]=CELLCENTERED)"
    write(31,*) hxx(1:nx+1,1:nz+1)
    write(31,*) hyy(1:nx+1,1:nz+1)
    write(31,*) hzz(1:nx+1,1:nz+1)
    write(31,*) huu(1:nx,1:nz)
    write(31,*) hwwma(1:nx,1:nz)
    write(31,*) hpp(1:nx,1:nz)
    write(31,*) hss(1:nx,1:nz)
    write(31,*) hmmiu(1:nx,1:nz)
    close(31)
    deallocate(hx)
    deallocate(hr)
    deallocate(sss)
    deallocate(hs)
    deallocate(hmiu)
    deallocate(hu)
    deallocate(hp)
    deallocate(hwma)
    deallocate(huu)
    deallocate(hwwma)
    deallocate(hpp)
    deallocate(hss)
    deallocate(hmmiu)
    end subroutine span
      
subroutine chord
    use global
    implicit none
    real(8) :: temp,tem,x1,x2,x3,r1,t1,t2,qqw,qqr
    real(8),allocatable :: huu(:,:),hpp(:,:),hss(:,:),hwwma(:,:),hmmiu(:,:)
    
    allocate(huu(1:ny,1:nz))
    allocate(hpp(1:ny,1:nz))
    allocate(hss(1:ny,1:nz))
    allocate(hwwma(1:ny,1:nz))
    allocate(hmmiu(1:ny,1:nz))
    do k=1,nz
        do j=1,ny
            y1=yy0(i,j,k)
            z1=zz0(i,j,k)
            huu(j,k)=Q12(i,j,k)/Q11(i,j,k)
            vy=Q13(i,j,k)/Q11(i,j,k)
            vz=Q14(i,j,k)/Q11(i,j,k)
            qq2=huu(j,k)*huu(j,k)+vy*vy+vz*vz!绝对速度平方
            wx=huu(j,k)
            wy=vy+rpm(rm)*z1
            wz=vz-rpm(rm)*y1
            qqw=wx*wx+wy*wy+wz*wz          !相对速度平方
            hpp(j,k)=0.4d0*(Q15(i,j,k)-0.5d0*Q11(i,j,k)*qq2)!静压
            hss(j,k)=ccv*log((hpp(j,k)/ppin)/((Q11(i,j,k)/dimin)**1.4))
            a=1.4d0*hpp(j,k)/Q11(i,j,k)
            hwwma(j,k)=sqrt(qqw/a)  !相对马赫数
            tem=hpp(j,k)/(Q11(i,j,k)*rg)!静温
            cvl=cvl0*((tem/t0)**1.5)*(t0+ts)/(tem+ts)!分子粘度
            hmmiu(j,k)=Q16(i,j,k)/cvl           !涡粘度/分子粘度
        end do
    end do
    open(31,file=trim(adjustl(secname))//'-'//trim(adjustl(id_mm))//'period-'//trim(adjustl(id_m))//'.dat')
    write(31,*) "VARIABLES=x,y,z,pvx,wma,pressure,entropy,viscosity"  !抬头
    write(31,"(1x,A,i4,2x,A,i4,2x,A,i4,2x,A)") "ZONE I=",1,"J=",ny+1,"K=",nz+1,"DATAPACKING=BLOCK, VARLOCATION=([4-8]=CELLCENTERED)"
    write(31,*) x(i,1:ny+1,1:nz+1)!用初始时刻的坐标值，相当于相对坐标系，所以失速团逆着旋转方向旋转
    write(31,*) y(i,1:ny+1,1:nz+1)
    write(31,*) z(i,1:ny+1,1:nz+1)
    write(31,*) huu(1:ny,1:nz)
    write(31,*) hwwma(1:ny,1:nz)
    write(31,*) hpp(1:ny,1:nz)
    write(31,*) hss(1:ny,1:nz)
    write(31,*) hmmiu(1:ny,1:nz)
    close(31)
    deallocate(huu)
    deallocate(hwwma)
    deallocate(hpp)
    deallocate(hss)
    deallocate(hmmiu)
    end subroutine chord
    
subroutine wl(medx,medr,spax,spar,n1,n2)
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

