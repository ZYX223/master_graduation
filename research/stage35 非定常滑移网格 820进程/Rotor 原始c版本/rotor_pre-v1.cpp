#include"rotor.h"
#include "mpi.h"
#include<cmath>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include <iomanip>
#include<algorithm>
//#include <unistd.h>
#include <iomanip>


namespace cfd
{
	int nb, nt, ng, nxm1, nym1, nzm1, nng0, nitt, nng, ign, it00, lbb, nxm, nym, nzm, ibm, itm, jbm, jtm;
	int layer;
	int* cg;
	int nx, ny, nz, ib, it, jb, jt, n, ci, cycle_by;
	double vxx, vrr, vtt, vee, y1, z1, rr, vx, vy, vz, wx, wy, wz, ama, dim, en, pp, rpm, ma;
	double c2, cfl, a4, a2, beta1, beta2;
	double time_begin, time_end, time_t;
	double cor, sir;
	double pi = 3.141592654;
	int* nnx, *nny, *nnz, *nib, *nit, *njb, *njt;
	double ***q01, ***q02, ***q03, ***q04, ***q05, ***q06;
	double ****q31, ****q32, ****q33, ****q34, ****q35, ****q36;
	double ****q11, ****q12, ****q13, ****q14, ****q15, ****q16;
	double ***av1, ***av2, ***av3, ***av4, ***av5, ***av6;
	double ***qc1, ***qc2, ***qc3, ***qc4, ***qc5, ***qc6;
	double ***qv1, ***qv2, ***qv3, ***qv4, ***qv5, ***qv6;
	double**** gradfi, ****gradfj, ****gradfk, ****gradc, ****gradcs;
	double*** ts1, ***ts2, ***ts3, ***ts4, ***ts5, ***ts6;
	double*** py1, ***py2, ***py3, ***py4, ***py5, ***py6;
	double*** rr1, ***rr2, ***rr3, ***rr4, ***rr5;
	double*** qp1, ***qp2, ***qp3, ***qp4, ***qp5;
	double ***xf, ***yf, ***zf, ***x, ***y, ***z, ***xx00, ***yy00, ***zz00;
	double ****xx1, ****yy1, ****zz1, ****xx2, ****yy2, ****zz2, ****xx3, ****yy3, ****zz3, ****xx, ****yy, ****zz;
	//double*** xx01, ***yy01, ***zz01, ***xx02, ***yy02, ***zz02, ***xx03, ***yy03, ***zz03;
	double ***xx0, ***yy0, ***zz0;
	double**** s1xn, ****s1yn, ****s1zn, ****s2xn, ****s2yn, ****s2zn, ****s3xn, ****s3yn, ****s3zn, ****vvn;
	//double ***s1x, ***s1y, ***s1z, ***s2x, ***s2y, ***s2z, ***s3x, ***s3y, ***s3z;
	//double ***vv;
	double ***pvx, ***pvy, ***pvz, ***vth, ***vre, ***p, ***t, ***time, ***wma;
	double ***dmini;
	double ***sri, ***srj, ***srk;
	double** dm;
	double* rms;
	double**** betaxn, ****betayn, ****betazn, ****petn, ****pebn, ****turin, ****hatn;
	double **betax, **betay, **betaz, **pet, **peb, **turi, **hat;
	double ta, timl, pt, ht, rout, pb0, pb1, period, rmsm, rmsm0, rmsmmax;
	double cvl0, t0, ts, cp, prt, prl, rg, cv1, cv2, kap, sigmav, cb1, cb2, cw1, cw2, cw3, cr1, cr2, cr3;
	string id_m;
	fstream f, f1, f2;
	int myid, myidl, myidr, numprocs, ierr, rc;
	MPI_Status status;

	std::string& trim(std::string &s)
	{
		if (s.empty())
		{
			return s;
		}
		s.erase(0, s.find_first_not_of(" "));
		s.erase(s.find_last_not_of(" ") + 1);
		return s;
	}
	//子例行程序

	//ini  读入控制参数和网格
	void ini()
	{
		double temp, cor1, sir1, t1, t2;

		//湍流计算所需参数
		cvl0 = 1.7161e-5;
		t0 = 273.16;
		ts = 110.4;
		rg = 287.0;
		cp = rg*1.4 / 0.4;

		prl = 0.72;
		prt = 0.9;

		cv1 = 7.1;
		cv2 = 5.0;

		kap = 0.41;

		cb1 = 0.1355;
		cb2 = 0.622;

		sigmav = 2.0 / 3.0;

		cw1 = cb1 / pow(kap, 2) + (1.0 + cb2) / sigmav;
		cw2 = 0.3;
		cw3 = 2.0;

		cr1 = 1.0;
		cr2 = 2.0;
		cr3 = 1.0;

		ifstream fin("ini3.dat"); //文件读取应检查文件是否打开成功 ，以及读写异常
		fin >> nt >> ng >> beta1 >> beta2 >> cfl >> a2 >> a4 >> ht >> pt >> pb1 >> c2 >> rmsm0;

		cg = (int*)malloc((ng + 1) * sizeof(int));

		for (int i = 1; i<ng + 1; i++)
		{
			double tmp;
			fin >> tmp;
			cg[i] = tmp;
		}

		fin >> nxm >> nym >> nzm >> ibm >> itm >> jbm >> jtm >> lbb >> rpm >> ma;
		rpm = rpm * pi / 30.0;
		fin.close();
		ifstream fin1("grid.dat");

		fin1 >> nxm1 >> nym1 >> nzm1;
		//初始化xf，yf，zf 行数
		xf = malloc3D<double>(nxm + 2, nym + 2, nzm + 2);
		yf = malloc3D<double>(nxm + 2, nym + 2, nzm + 2);
		zf = malloc3D<double>(nxm + 2, nym + 2, nzm + 2);



		//数组赋值

		for (int i = 0; i<nzm1 + 1; i++)
			for (int j = 0; j<nym1 + 1; j++)
				for (int k = 0; k<nxm1 + 1; k++)
				{
					if (k == 0 || j == 0 || i == 0)
					{
						xf[i][j][k] = yf[i][j][k] = zf[i][j][k] = 0;
					}
					else
					{
						fin1 >> xf[i][j][k] >> yf[i][j][k] >> zf[i][j][k];
					}
				}

		fin1.close();
		//由单叶道数据旋转得到整圈叶排数据
		nx = nxm;
		ny = nym;
		nz = nzm;

		temp = 2.0 / double(lbb)*pi*myid;

		cor1 = cos(temp);
		sir1 = sin(temp);
		printf("nt=%d,rmsm0=%f,ma=%f,nx=%d\n",nt,rmsm0,ma,nx);
		//		for (int k = 1; k <nz + 2; k++)

		//		for (int j = 1; j < ny + 2; j++)

		for (int i = 1; i< nz + 2; i++)
			for (int j = 1; j<ny + 2; j++)
				for (int k = 1; k<nx + 2; k++)
				{
					t1 = yf[i][j][k];
					t2 = zf[i][j][k];

					yf[i][j][k] = t1*cor1 - t2*sir1;
					zf[i][j][k] = t1*sir1 + t2*cor1;
				}

		ostringstream os;
		os << myid;
		id_m = os.str();
	}
	//allocation  申请内存

	void allocation()
	{
		nnx = (int *)malloc((ng + 1) * sizeof(int));
		nny = (int *)malloc((ng + 1) * sizeof(int));
		nnz = (int *)malloc((ng + 1) * sizeof(int));
		nib = (int *)malloc((ng + 1) * sizeof(int));
		nit = (int *)malloc((ng + 1) * sizeof(int));
		njb = (int *)malloc((ng + 1) * sizeof(int));
		njt = (int *)malloc((ng + 1) * sizeof(int));

		dm = malloc2D<double>(nt + 1, nt + 1);

		q11 = malloc4D<double>(nx + 2, ny + 2, nz + 2, nt + 1);
		q12 = malloc4D<double>(nx + 2, ny + 2, nz + 2, nt + 1);
		q13 = malloc4D<double>(nx + 2, ny + 2, nz + 2, nt + 1);
		q14 = malloc4D<double>(nx + 2, ny + 2, nz + 2, nt + 1);
		q15 = malloc4D<double>(nx + 2, ny + 2, nz + 2, nt + 1);
		q16 = malloc4D<double>(nx + 2, ny + 2, nz + 2, nt + 1);

		q31 = malloc4D<double>(nx + 1, ny + 1, nz + 1, nt + 1);
		q32 = malloc4D<double>(nx + 1, ny + 1, nz + 1, nt + 1);
		q33 = malloc4D<double>(nx + 1, ny + 1, nz + 1, nt + 1);
		q34 = malloc4D<double>(nx + 1, ny + 1, nz + 1, nt + 1);
		q35 = malloc4D<double>(nx + 1, ny + 1, nz + 1, nt + 1);
		q36 = malloc4D<double>(nx + 1, ny + 1, nz + 1, nt + 1);

		x = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		y = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		z = malloc3D<double>(nx + 2, ny + 2, nz + 2);

		xx00 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		yy00 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		zz00 = malloc3D<double>(nx + 1, ny + 1, nz + 1);

		dmini = malloc3D<double>(nx + 1, ny + 1, nz + 1);

		betaxn = malloc4D<double>(ny + 1, nz + 1, ng + 1, nt + 1);
		betayn = malloc4D<double>(ny + 1, nz + 1, ng + 1, nt + 1);
		betazn = malloc4D<double>(ny + 1, nz + 1, ng + 1, nt + 1);
		hatn = malloc4D<double>(ny + 1, nz + 1, ng + 1, nt + 1);
		petn = malloc4D<double>(ny + 1, nz + 1, ng + 1, nt + 1);
		turin = malloc4D<double>(ny + 1, nz + 1, ng + 1, nt + 1);
		pebn = malloc4D<double>(ny + 1, nz + 1, ng + 1, nt + 1);

		betax = malloc2D<double>(ny + 1, nz + 1);
		betay = malloc2D<double>(ny + 1, nz + 1);
		betaz = malloc2D<double>(ny + 1, nz + 1);
		hat = malloc2D<double>(ny + 1, nz + 1);
		pet = malloc2D<double>(ny + 1, nz + 1);
		turi = malloc2D<double>(ny + 1, nz + 1);
		peb = malloc2D<double>(ny + 1, nz + 1);

		ts1 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		ts2 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		ts3 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		ts4 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		ts5 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		ts6 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		//vv = malloc3D<double>(nx + 1, ny + 1, nz + 1);

		q01 = malloc3D<double>(nx + 1, ny + 1, nz + 1);

		q02 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		q03 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		q04 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		q05 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		q06 = malloc3D<double>(nx + 1, ny + 1, nz + 1);

		py1 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		py2 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		py3 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		py4 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		py5 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		py6 = malloc3D<double>(nx + 1, ny + 1, nz + 1);

		rr1 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		rr2 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		rr3 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		rr4 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		rr5 = malloc3D<double>(nx + 1, ny + 1, nz + 1);

		qp1 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		qp2 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		qp3 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		qp4 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		qp5 = malloc3D<double>(nx + 1, ny + 1, nz + 1);

		time = malloc3D<double>(nx + 1, ny + 1, nz + 1);

		s2xn = malloc4D<double>(nx + 3, ny + 1, nz + 1, ng + 1);
		s2yn = malloc4D<double>(nx + 3, ny + 1, nz + 1, ng + 1);
		s2zn = malloc4D<double>(nx + 3, ny + 1, nz + 1, ng + 1);
		xx2 = malloc4D<double>(nx + 3, ny + 1, nz + 1, ng + 1);

		yy2 = malloc4D<double>(nx + 2, ny + 1, nz + 1, ng + 1);
		zz2 = malloc4D<double>(nx + 2, ny + 1, nz + 1, ng + 1);

		s3xn = malloc4D<double>(nx + 1, ny + 3, nz + 1, ng + 1);
		s3yn = malloc4D<double>(nx + 1, ny + 3, nz + 1, ng + 1);
		s3zn = malloc4D<double>(nx + 1, ny + 3, nz + 1, ng + 1);

		xx3 = malloc4D<double>(nx + 1, ny + 3, nz + 1, ng + 1);
		yy3 = malloc4D<double>(nx + 1, ny + 2, nz + 1, ng + 1);
		zz3 = malloc4D<double>(nx + 1, ny + 2, nz + 1, ng + 1);


		s1xn = malloc4D<double>(nx + 1, ny + 1, nz + 3, ng + 1);
		s1yn = malloc4D<double>(nx + 1, ny + 1, nz + 3, ng + 1);
		s1zn = malloc4D<double>(nx + 1, ny + 1, nz + 3, ng + 1);

		xx1 = malloc4D<double>(nx + 1, ny + 1, nz + 2, ng + 1);
		yy1 = malloc4D<double>(nx + 1, ny + 1, nz + 2, ng + 1);
		zz1 = malloc4D<double>(nx + 1, ny + 1, nz + 2, ng + 1);

		vvn = malloc4D<double>(nx + 1, ny + 1, nz + 1, ng + 1);
		xx = malloc4D<double>(nx + 1, ny + 1, nz + 1, ng + 1);

		yy = malloc4D<double>(nx + 1, ny + 1, nz + 2, ng + 1);
		zz = malloc4D<double>(nx + 1, ny + 1, nz + 2, ng + 1);

		/*s1x = malloc3D<double>(nx + 1, ny + 1, nz + 3);
		s1y = malloc3D<double>(nx + 1, ny + 1, nz + 3);
		s1z = malloc3D<double>(nx + 1, ny + 1, nz + 3);

		s2x = malloc3D<double>(nx + 3, ny + 1, nz + 1);
		s2y = malloc3D<double>(nx + 3, ny + 1, nz + 1);
		s2z = malloc3D<double>(nx + 3, ny + 1, nz + 1);

		s3x = malloc3D<double>(nx + 1, ny + 3, nz + 1);
		s3y = malloc3D<double>(nx + 1, ny + 3, nz + 1);
		s3z= malloc3D<double>(nx + 1, ny + 3, nz + 1);
		//*/
		/*xx02 = malloc3D<double>(nx + 2, ny + 1, nz + 1);
		yy02 = malloc3D<double>(nx + 2, ny + 1, nz + 1);
		zz02 = malloc3D<double>(nx + 2, ny + 1, nz + 1);

		xx01 = malloc3D<double>(nx + 1, ny + 1, nz + 2);
		yy01 = malloc3D<double>(nx + 1, ny + 1, nz + 2);
		zz01 = malloc3D<double>(nx + 1, ny + 1, nz + 2);

		xx03 = malloc3D<double>(nx + 1, ny + 2, nz + 1);
		yy03 = malloc3D<double>(nx + 1, ny + 2, nz + 1);
		zz03 = malloc3D<double>(nx + 1, ny + 2, nz + 1);
		//*/
		xx0 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		yy0 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		zz0 = malloc3D<double>(nx + 2, ny + 2, nz + 2);

		av1 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		av2 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		av3 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		av4 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		av5 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		av6 = malloc3D<double>(nx + 2, ny + 2, nz + 2);

		qc1 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		qc2 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		qc3 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		qc4 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		qc5 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		qc6 = malloc3D<double>(nx + 2, ny + 2, nz + 2);

		qv1 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		qv2 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		qv3 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		qv4 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		qv5 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		qv6 = malloc3D<double>(nx + 2, ny + 2, nz + 2);

		pvx = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		pvy = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		pvz = malloc3D<double>(nx + 2, ny + 2, nz + 2);

		vth = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		vre = malloc3D<double>(nx + 2, ny + 2, nz + 2);

		p = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		t = malloc3D<double>(nx + 2, ny + 2, nz + 2);

		wma = malloc3D<double>(nx + 2, ny + 2, nz + 2);

		sri = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		srj = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		srk = malloc3D<double>(nx + 2, ny + 2, nz + 2);

		gradfi = malloc4D<double>(nx + 2, ny + 2, nz + 2, 16);

		gradfj = malloc4D<double>(nx + 2, ny + 2, nz + 2, 16);

		gradfk = malloc4D<double>(nx + 2, ny + 2, nz + 2, 16);

		gradc = malloc4D<double>(nx + 2, ny + 2, nz + 2, 13);

		gradcs = malloc4D<double>(nx + 2, ny + 2, nz + 2, 10);

		rms = (double*)malloc((nt + 1) * sizeof(double));
	}

	void freealloc()
	{
		free(nnx);
		free(nny);
		free(nnz);
		free(nib);
		free(nit);
		free(njb);
		free(njt);
		free2D<double>(dm);
		free4D<double>(q11);
		free4D<double>(q12);
		free4D<double>(q13);
		free4D<double>(q14);
		free4D<double>(q15);
		free4D<double>(q16);
		free4D<double>(q31);
		free4D<double>(q32);
		free4D<double>(q33);
		free4D<double>(q34);
		free4D<double>(q35);
		free4D<double>(q36);
		free3D<double>(x);
		free3D<double>(y);
		free3D<double>(z);
		free3D<double>(xx00);
		free3D<double>(yy00);
		free3D<double>(zz00);
		free3D<double>(dmini);
		free4D<double>(betaxn);
		free4D<double>(betayn);
		free4D<double>(betazn);
		free4D<double>(hatn);
		free4D<double>(petn);
		free4D<double>(turin);
		free4D<double>(pebn);
		free2D<double>(betax);
		free2D<double>(betay);
		free2D<double>(betaz);
		free2D<double>(hat);
		free2D<double>(pet);
		free2D<double>(turi);
		free2D<double>(peb);
		free3D<double>(ts1);
		free3D<double>(ts2);
		free3D<double>(ts3);
		free3D<double>(ts4);
		free3D<double>(ts5);
		free3D<double>(ts6);
		free3D<double>(q01);
		free3D<double>(q02);
		free3D<double>(q03);
		free3D<double>(q04);
		free3D<double>(q05);
		free3D<double>(q06);
		free3D<double>(py1);
		free3D<double>(py2);
		free3D<double>(py3);
		free3D<double>(py4);
		free3D<double>(py5);
		free3D<double>(py6);
		free3D<double>(qp1);
		free3D<double>(qp2);
		free3D<double>(qp3);
		free3D<double>(qp4);
		free3D<double>(qp5);
		free3D<double>(rr1);
		free3D<double>(rr2);
		free3D<double>(rr3);
		free3D<double>(rr4);
		free3D<double>(rr5);
		free3D<double>(time);
		free4D<double>(s2xn);
		free4D<double>(s2yn);
		free4D<double>(s2zn);
		free4D<double>(xx2);
		free4D<double>(yy2);
		free4D<double>(zz2);
		free4D<double>(s3xn);
		free4D<double>(s3yn);
		free4D<double>(s3zn);
		free4D<double>(xx3);
		free4D<double>(yy3);
		free4D<double>(zz3);
		free4D<double>(s1xn);
		free4D<double>(s1yn);
		free4D<double>(s1zn);
		free4D<double>(xx1);
		free4D<double>(yy1);
		free4D<double>(zz1);
		free4D<double>(xx);
		free4D<double>(yy);
		free4D<double>(zz);
		free4D<double>(vvn);
		free3D<double>(xx0);
		free3D<double>(yy0);
		free3D<double>(zz0);
		free3D<double>(qc1);
		free3D<double>(qc2);
		free3D<double>(qc3);
		free3D<double>(qc4);
		free3D<double>(qc5);
		free3D<double>(qc6);
		free3D<double>(qv1);
		free3D<double>(qv2);
		free3D<double>(qv3);
		free3D<double>(qv4);
		free3D<double>(qv5);
		free3D<double>(qv6);
		free3D<double>(pvx);
		free3D<double>(pvy);
		free3D<double>(pvz);
		free3D<double>(vth);
		free3D<double>(vre);
		free3D<double>(p);
		free3D<double>(t);
		free3D<double>(wma);
		free3D<double>(sri);
		free3D<double>(srj);
		free3D<double>(srk);
		free4D<double>(gradfi);
		free4D<double>(gradfj);
		free4D<double>(gradfk);
		free4D<double>(gradc);
		free4D<double>(gradcs);
		free(rms);
	}

	// shu  多重网格数

	void shu()
	{
		double kk;

		itm = itm - 1;
		jtm = jtm - 1;

		for (int j = 1; j < ng + 1; j++)
		{
			kk = pow(2.0, (ng - j));

			nnx[j] = nx / kk;
			nny[j] = ny / kk;
			nnz[j] = nz / kk;

			nib[j] = (ibm - 1) / kk + 1;
			nit[j] = itm / kk;
			njb[j] = (jbm - 1) / kk + 1;
			njt[j] = jtm / kk;
		}
	}


	void ini2()
	{
		double a;
		//double sir,cor;
		int kk, ifine, jfine, kfine;

		//初始化 nng0，表示粗网格

		nng0 = 1;

		nx = nnx[nng0];
		ny = nny[nng0];
		nz = nnz[nng0];

		kk = pow(2.0, (ng - nng0));

		/*	for (int i = 1; i<nx + 2; i++)

		{

		ifine = kk*(i - 1) + 1;

		for (int j = 1; j<ny + 2; j++)

		{

		jfine = kk*(j - 1) + 1;

		for (int k = 1; i<nz + 2; k++)

		{

		kfine = kk*(k - 1) + 1;

		x[i][j][k] = xf[ifine][jfine][kfine];

		y[i][j][k] = yf[ifine][jfine][kfine];

		z[i][j][k] = zf[ifine][jfine][kfine];

		}

		}

		}*/

		//输出全0
		for (int k = 1; k<nz + 2; k++)
		{
			kfine = kk*(k - 1) + 1;
			for (int j = 1; j<ny + 2; j++)
			{
				jfine = kk*(j - 1) + 1;
				for (int i = 1; i<nx + 2; i++)
				{
					ifine = kk*(i - 1) + 1;
					x[k][j][i] = xf[kfine][jfine][ifine];
					y[k][j][i] = yf[kfine][jfine][ifine];
					z[k][j][i] = zf[kfine][jfine][ifine];
				}
			}
		}
		vxx = 1.0;
		vtt = vxx*beta1;
		vrr = vxx*beta2;
		vee = sqrt(pow(vxx, 2) + pow(vrr, 2) + pow(vtt, 2));

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
			{
				y1 = 0.25*(y[k][j][1] + y[k][j + 1][1] + y[k + 1][j][1] + y[k + 1][j + 1][1]);
				z1 = 0.25*(z[k][j][1] + z[k][j + 1][1] + z[k + 1][j][1] + z[k + 1][j + 1][1]);

				rr = sqrt(y1*y1 + z1*z1);

				sir = z1 / rr;
				cor = y1 / rr;

				for (int i = 1; i<nt + 1; i++)
				{
					betaxn[i][nng0][k][j] = vxx / vee;
					betayn[i][nng0][k][j] = (vrr*cor - vtt*sir) / vee;
					betazn[i][nng0][k][j] = (vrr*sir + vtt*cor) / vee;
				}
			}
		for (int k = 1; k<nt + 1; k++)
			for (int j = 1; j<nz + 1; j++)
				for (int i = 1; i<ny + 1; i++)
				{
					hatn[k][nng0][j][i] = ht;
					petn[k][nng0][j][i] = pt;
				}

		//初始化总迭代次数

		nitt = 0;
		for (int n = 1; n<nt + 1; n++)
			for (int k = 1; k<nz + 1; k++)
				for (int j = 1; j<ny + 1; j++)
					for (int i = 1; i<nx + 1; i++)
					{
						rout = 3.5*petn[n][nng0][k][j] / hatn[n][nng0][k][j];
						pp = petn[n][nng0][k][j] * pow((1.0 + 0.2*ma*ma), -3.5);
						dim = rout*pow((1.0 + 0.2*ma*ma), -2.5);
						a = sqrt(1.4*pp / dim);
						en = 2.5*pp + 0.5*dim*ma*ma*a*a;

						q11[n][k][j][i] = dim;
						q12[n][k][j][i] = dim*a*ma*betaxn[n][nng0][k][j];
						q13[n][k][j][i] = dim*a*ma*betayn[n][nng0][k][j];
						q14[n][k][j][i] = dim*a*ma*betazn[n][nng0][k][j];
						q15[n][k][j][i] = en;
						q16[n][k][j][i] = 200 * cvl0;
					}

		string file = "error-" + trim(id_m) + "myid.dat";
		f.open(file.c_str(),std::fstream::out | std::fstream::app);

		if (myid == 0)
		{
			f1.open("convergence-inflow.dat", std::fstream::out | std::fstream::app);
			f2.open("convergence-outflow.dat", std::fstream::out | std::fstream::app);
		}
	}

	//时间谱系数

	void dodd()
	{
		period = 2.0*pi / rpm;

		for (int i = 1; i<nt + 1; i++)
			for (int j = 1; j<nt + 1; j++)
			{
				dm[i][j] = 0;
			}

		if (nt % 2 == 1)
		{
			for (int i = 1; i<nt + 1; i++)
				for (int j = 1; j<nt + 1; j++)
				{
					if (i != j)
						dm[i][j] = pi / period*pow(-1.0, (i - j)) / sin(pi*double((i - j)) / double(nt));
					//	dm[i][j] = pi / period*pow(-1.0, (j - i)) / sin(pi*double((j - i)) / double(nt));
				}
		}
		else
		{
			for (int i = 1; i<nt + 1; i++)
				for (int j = 1; j<nt + 1; j++)
				{
					if (i != j)
						dm[i][j] = pi / period*pow(-1.0, (i - j)) / tan(pi*double((i - j)) / double(nt));
					//	dm[i][j] = pi / period*pow(-1.0, (j - i)) / tan(pi*double((j - i)) / double(nt));
				}
		}
	}

	void distance()
	{
		int ii, jj, kk, iil, iir, kkl, kkr, iib, iit, jjb, jjt;
		double d, dy1, dy2, dy, dz1, dz2, dz;

		nx = nnx[ng];
		ny = nny[ng];
		nz = nnz[ng];

		ib = nib[ng];
		it = nit[ng];
		jb = njb[ng];
		jt = njt[ng];

		vector<vector<double> > xwu(nz + 1, vector<double>(nx + 1)), ywu(nz + 1, vector<double>(nx + 1)), zwu(nz + 1, vector<double>(nx + 1));
		vector<vector<double> > xwd(nz + 1, vector<double>(nx + 1)), ywd(nz + 1, vector<double>(nx + 1)), zwd(nz + 1, vector<double>(nx + 1));
		vector<vector<double> > xwf(jt + 1, vector<double>(it + 1)), ywf(jt + 1, vector<double>(it + 1)), zwf(jt + 1, vector<double>(it + 1));
		vector<vector<double> > xwb(jt + 1, vector<double>(it + 1)), ywb(jt + 1, vector<double>(it + 1)), zwb(jt + 1, vector<double>(it + 1));

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					xx00[k][j][i] = 0.125*(xf[k][j][i] + xf[k][j][i + 1] + xf[k + 1][j][i] + xf[k + 1][j][i + 1]
						+ xf[k][j + 1][i] + xf[k][j + 1][i + 1] + xf[k + 1][j + 1][i] + xf[k + 1][j + 1][i + 1]);

					yy00[k][j][i] = 0.125*(yf[k][j][i] + yf[k][j][i + 1] + yf[k + 1][j][i] + yf[k + 1][j][i + 1]
						+ yf[k][j + 1][i] + yf[k][j + 1][i + 1] + yf[k + 1][j + 1][i] + yf[k + 1][j + 1][i + 1]);

					zz00[k][j][i] = 0.125*(zf[k][j][i] + zf[k][j][i + 1] + zf[k + 1][j][i] + zf[k + 1][j][i + 1]
						+ zf[k][j + 1][i] + zf[k][j + 1][i + 1] + zf[k + 1][j + 1][i] + zf[k + 1][j + 1][i + 1]);
				}


		for (int k = 1; k<nz + 1; k++)
			for (int i = 1; i<nx + 1; i++)
			{
				xwd[k][i] = 0.25*(xf[k][1][i] + xf[k][1][i + 1] + xf[k + 1][1][i] + xf[k + 1][1][i + 1]);
				ywd[k][i] = 0.25*(yf[k][1][i] + yf[k][1][i + 1] + yf[k + 1][1][i] + yf[k + 1][1][i + 1]);
				zwd[k][i] = 0.25*(zf[k][1][i] + zf[k][1][i + 1] + zf[k + 1][1][i] + zf[k + 1][1][i + 1]);



				xwu[k][i] = 0.25*(xf[k][ny + 1][i] + xf[k][ny + 1][i + 1] + xf[k + 1][ny + 1][i] + xf[k + 1][ny + 1][i + 1]);
				ywu[k][i] = 0.25*(yf[k][ny + 1][i] + yf[k][ny + 1][i + 1] + yf[k + 1][ny + 1][i] + yf[k + 1][ny + 1][i + 1]);
				zwu[k][i] = 0.25*(zf[k][ny + 1][i] + zf[k][ny + 1][i + 1] + zf[k + 1][ny + 1][i] + zf[k + 1][ny + 1][i + 1]);
			}

		for (int j = jb; j<jt + 1; j++)
			for (int i = ib; i<it + 1; i++)
			{
				xwf[j][i] = 0.25*(xf[1][j][i] + xf[1][j][i + 1] + xf[1][j + 1][i] + xf[1][j + 1][i + 1]);
				ywf[j][i] = 0.25*(yf[1][j][i] + yf[1][j][i + 1] + yf[1][j + 1][i] + yf[1][j + 1][i + 1]);
				zwf[j][i] = 0.25*(zf[1][j][i] + zf[1][j][i + 1] + zf[1][j + 1][i] + zf[1][j + 1][i + 1]);



				xwb[j][i] = 0.25*(xf[nz + 1][j][i] + xf[nz + 1][j][i + 1] + xf[nz + 1][j + 1][i] + xf[nz + 1][j + 1][i + 1]);
				ywb[j][i] = 0.25*(yf[nz + 1][j][i] + yf[nz + 1][j][i + 1] + yf[nz + 1][j + 1][i] + yf[nz + 1][j + 1][i + 1]);
				zwb[j][i] = 0.25*(zf[nz + 1][j][i] + zf[nz + 1][j][i + 1] + zf[nz + 1][j + 1][i] + zf[nz + 1][j + 1][i + 1]);
			}

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					dy = pow(10.0, 20);

					iil = i - 10;
					iil = max(iil, 1);
					iir = iil + 20;
					iir = min(iir, nx);



					kkl = k - 10;
					kkl = max(kkl, 1);
					kkr = kkl + 20;
					kkr = min(kkr, nz);



					for (kk = kkl; kk<kkr + 1; kk++)
						for (ii = iil; ii<iir + 1; ii++)
						{
							dy1 = pow((xx00[k][j][i] - xwd[kk][ii]), 2) + pow((yy00[k][j][i] - ywd[kk][ii]), 2) + pow((zz00[k][j][i] - zwd[kk][ii]), 2);
							dy2 = pow((xx00[k][j][i] - xwu[kk][ii]), 2) + pow((yy00[k][j][i] - ywu[kk][ii]), 2) + pow((zz00[k][j][i] - zwu[kk][ii]), 2);

							d = min(dy1, dy2);
							if (d<dy)
								dy = d;
						}



					dz = pow(10.0, 20);

					if (j<jb)
					{
						jjb = jb;
						jjt = jb + 20;
					}
					else if ((j >= jb) && (j <= jt))
					{
						jjb = j - 10;
						jjb = max(jjb, jb);
						jjt = jjb + 20;
						jjt = min(jjt, jt);
					}
					else
					{
						jjb = jt - 20;
						jjt = jt;
					}


					if (i<ib)
					{
						iib = ib;
						iit = ib + 20;
					}
					else if ((i >= ib) && (i <= it))
					{
						iib = i - 10;
						iib = max(iib, ib);
						iit = iib + 20;
						iit = min(iit, it);
					}
					else
					{
						iib = it - 20;
						iit = it;
					}


					for (jj = jjb; jj<jjt + 1; jj++)
						for (ii = iib; ii<iit + 1; ii++)
						{
							dz1 = (pow((xx00[k][j][i] - xwf[jj][ii]), 2) + pow((yy00[k][j][i] - ywf[jj][ii]), 2) + pow((zz00[k][j][i] - zwf[jj][ii]), 2));
							dz2 = (pow((xx00[k][j][i] - xwb[jj][ii]), 2) + pow((yy00[k][j][i] - ywb[jj][ii]), 2) + pow((zz00[k][j][i] - zwb[jj][ii]), 2));

							d = min(dz1, dz2);
							if (d<dz)
								dz = d;
						}
					dmini[k][j][i] = sqrt(min(dy, dz));

				}
	}

	void geo1(double*** x0, double*** y0, double*** z0,
		double**** s01x, double**** s01y, double**** s01z,
		double**** xx001, double**** yy001, double**** zz001,
		double**** s02x, double**** s02y, double**** s02z,
		double**** xx002, double**** yy002, double**** zz002,
		double**** s03x, double**** s03y, double**** s03z,
		double**** xx003, double**** yy003, double**** zz003,
		double**** vv0, double**** xx000, double**** yy000,
		double**** zz000)
	{

		double x24, y24, z24, x31, y31, z31;

		for (int k = 1; k<nz + 2; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					x24 = x0[k][j + 1][i] - x0[k][j][i + 1];
					x31 = x0[k][j + 1][i + 1] - x0[k][j][i];

					y24 = y0[k][j + 1][i] - y0[k][j][i + 1];
					y31 = y0[k][j + 1][i + 1] - y0[k][j][i];

					z24 = z0[k][j + 1][i] - z0[k][j][i + 1];
					z31 = z0[k][j + 1][i + 1] - z0[k][j][i];



					s01x[nng][k][j][i] = 0.5*(y24*z31 - z24*y31);
					s01y[nng][k][j][i] = 0.5*(z24*x31 - x24*z31);
					s01z[nng][k][j][i] = 0.5*(x24*y31 - y24*x31);



					xx001[nng][k][j][i] = 0.25*(x0[k][j][i] + x0[k][j + 1][i] + x0[k][j][i + 1] + x0[k][j + 1][i + 1]);
					yy001[nng][k][j][i] = 0.25*(y0[k][j][i] + y0[k][j + 1][i] + y0[k][j][i + 1] + y0[k][j + 1][i + 1]);
					zz001[nng][k][j][i] = 0.25*(z0[k][j][i] + z0[k][j + 1][i] + z0[k][j][i + 1] + z0[k][j + 1][i + 1]);

				}

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 2; i++)
				{
					x24 = x0[k + 1][j + 1][i] - x0[k][j][i];
					x31 = x0[k][j + 1][i] - x0[k + 1][j][i];

					y24 = y0[k + 1][j + 1][i] - y0[k][j][i];
					y31 = y0[k][j + 1][i] - y0[k + 1][j][i];

					z24 = z0[k + 1][j + 1][i] - z0[k][j][i];
					z31 = z0[k][j + 1][i] - z0[k + 1][j][i];



					s02x[nng][k][j][i] = 0.5*(y24*z31 - z24*y31);
					s02y[nng][k][j][i] = 0.5*(z24*x31 - x24*z31);
					s02z[nng][k][j][i] = 0.5*(x24*y31 - y24*x31);



					xx002[nng][k][j][i] = 0.25*(x0[k][j][i] + x0[k][j + 1][i] + x0[k + 1][j][i] + x0[k + 1][j + 1][i]);
					yy002[nng][k][j][i] = 0.25*(y0[k][j][i] + y0[k][j + 1][i] + y0[k + 1][j][i] + y0[k + 1][j + 1][i]);
					zz002[nng][k][j][i] = 0.25*(z0[k][j][i] + z0[k][j + 1][i] + z0[k + 1][j][i] + z0[k + 1][j + 1][i]);

				}

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 2; j++)
				for (int i = 1; i<nx + 1; i++)
				{

					x24 = x0[k + 1][j][i + 1] - x0[k][j][i];
					x31 = x0[k + 1][j][i] - x0[k][j][i + 1];

					y24 = y0[k + 1][j][i + 1] - y0[k][j][i];
					y31 = y0[k + 1][j][i] - y0[k][j][i + 1];

					z24 = z0[k + 1][j][i + 1] - z0[k][j][i];
					z31 = z0[k + 1][j][i] - z0[k][j][i + 1];



					s03x[nng][k][j][i] = 0.5*(y24*z31 - z24*y31);
					s03y[nng][k][j][i] = 0.5*(z24*x31 - x24*z31);
					s03z[nng][k][j][i] = 0.5*(x24*y31 - y24*x31);



					xx003[nng][k][j][i] = 0.25*(x0[k][j][i] + x0[k][j][i + 1] + x0[k + 1][j][i] + x0[k + 1][j][i + 1]);
					yy003[nng][k][j][i] = 0.25*(y0[k][j][i] + y0[k][j][i + 1] + y0[k + 1][j][i] + y0[k + 1][j][i + 1]);
					zz003[nng][k][j][i] = 0.25*(z0[k][j][i] + z0[k][j][i + 1] + z0[k + 1][j][i] + z0[k + 1][j][i + 1]);

				}

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{

					x24 = x0[k + 1][j + 1][i + 1] - x0[k][j][i];
					y24 = y0[k + 1][j + 1][i + 1] - y0[k][j][i];
					z24 = z0[k + 1][j + 1][i + 1] - z0[k][j][i];



					vv0[nng][k][j][i] = -(x24*(s01x[nng][k][j][i] + s02x[nng][k][j][i] + s03x[nng][k][j][i]) +
						y24*(s01y[nng][k][j][i] + s02y[nng][k][j][i] + s03y[nng][k][j][i]) +
						z24*(s01z[nng][k][j][i] + s02z[nng][k][j][i] + s03z[nng][k][j][i])) / 3.0;



					xx000[nng][k][j][i] = 0.125*(x0[k][j][i] + x0[k][j][i + 1] + x0[k + 1][j][i] + x0[k + 1][j][i + 1] +
						x0[k][j + 1][i] + x0[k][j + 1][i + 1] + x0[k + 1][j + 1][i] + x0[k + 1][j + 1][i + 1]);

					yy000[nng][k][j][i] = 0.125*(y0[k][j][i] + y0[k][j][i + 1] + y0[k + 1][j][i] + y0[k + 1][j][i + 1] +
						y0[k][j + 1][i] + y0[k][j + 1][i + 1] + y0[k + 1][j + 1][i] + y0[k + 1][j + 1][i + 1]);

					zz000[nng][k][j][i] = 0.125*(z0[k][j][i] + z0[k][j][i + 1] + z0[k + 1][j][i] + z0[k + 1][j][i + 1] +
						z0[k][j + 1][i] + z0[k][j + 1][i + 1] + z0[k + 1][j + 1][i] + z0[k + 1][j + 1][i + 1]);
				}

		for (int j = 1; j<ny + 1; j++)
			for (int i = 1; i<nx + 1; i++)
			{
				s01x[nng][0][j][i] = s01x[nng][1][j][i];
				s01y[nng][0][j][i] = s01y[nng][1][j][i];
				s01z[nng][0][j][i] = s01z[nng][1][j][i];

				s01x[nng][nz + 2][j][i] = s01x[nng][nz + 1][j][i];
				s01y[nng][nz + 2][j][i] = s01y[nng][nz + 1][j][i];
				s01z[nng][nz + 2][j][i] = s01z[nng][nz + 1][j][i];
			}
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
			{
				s02x[nng][k][j][0] = s02x[nng][k][j][1];
				s02y[nng][k][j][0] = s02y[nng][k][j][1];
				s02z[nng][k][j][0] = s02z[nng][k][j][1];

				s02x[nng][k][j][nx + 2] = s02x[nng][k][j][nx + 1];
				s02y[nng][k][j][nx + 2] = s02y[nng][k][j][nx + 1];
				s02z[nng][k][j][nx + 2] = s02z[nng][k][j][nx + 1];
			}

		//j-direction 边界值设定


		for (int k = 1; k<nz + 1; k++)
			for (int i = 1; i<nx + 1; i++)
			{
				s03x[nng][k][0][i] = s03x[nng][k][1][i];
				s03y[nng][k][0][i] = s03y[nng][k][1][i];
				s03z[nng][k][0][i] = s03z[nng][k][1][i];

				s03x[nng][k][ny + 2][i] = s03x[nng][k][ny + 1][i];
				s03y[nng][k][ny + 2][i] = s03y[nng][k][ny + 1][i];
				s03z[nng][k][ny + 2][i] = s03z[nng][k][ny + 1][i];
			}
	}

	void geo()
	{
		int ifine, jfine, kfine, kk;
		double val1, val2, val3, val4, val5, val6, val7, val8, val, temp, cor1, sir1;

		nx = nnx[nng];
		ny = nny[nng];
		nz = nnz[nng];
		kk = pow(2.0, ng - nng);
		for (int k = 1; k<nz + 2; k++)
		{
			kfine = kk*(k - 1) + 1;
			for (int j = 1; j<ny + 2; j++)
			{
				jfine = kk*(j - 1) + 1;
				for (int i = 1; i<nx + 2; i++)
				{
					ifine = kk*(i - 1) + 1;

					x[k][j][i] = xf[kfine][jfine][ifine];
					y[k][j][i] = yf[kfine][jfine][ifine];
					z[k][j][i] = zf[kfine][jfine][ifine];
				}
			}
		}
		//由顶点值计算面中心坐标值

		geo1(x, y, z, s1xn, s1yn, s1zn, xx1, yy1, zz1, s2xn, s2yn, s2zn, xx2, yy2, zz2, s3xn, s3yn, s3zn, xx3, yy3, zz3, vvn, xx, yy, zz);


		for (ign = nng - 1; ign >= 1; ign--)
		{
			nx = nnx[ign];
			ny = nny[ign];
			nz = nnz[ign];

			for (int k = 1; k<nz + 1; k++)
			{
				kfine = 2 * k - 1;
				for (int j = 1; j<ny + 1; j++)
				{
					jfine = 2 * j - 1;
					for (int i = 1; i<nx + 2; i++)
					{
						ifine = 2 * i - 1;

						val1 = sqrt(pow(s2xn[ign + 1][kfine][jfine][ifine], 2) + pow(s2yn[ign + 1][kfine][jfine][ifine], 2) + pow(s2zn[ign + 1][kfine][jfine][ifine], 2));
						val2 = sqrt(pow(s2xn[ign + 1][kfine][jfine + 1][ifine], 2) + pow(s2yn[ign + 1][kfine][jfine + 1][ifine], 2) + pow(s2zn[ign + 1][kfine][jfine + 1][ifine], 2));
						val3 = sqrt(pow(s2xn[ign + 1][kfine + 1][jfine][ifine], 2) + pow(s2yn[ign + 1][kfine + 1][jfine][ifine], 2) + pow(s2zn[ign + 1][kfine + 1][jfine][ifine], 2));
						val4 = sqrt(pow(s2xn[ign + 1][kfine + 1][jfine + 1][ifine], 2) + pow(s2yn[ign + 1][kfine + 1][jfine + 1][ifine], 2) + pow(s2zn[ign + 1][kfine + 1][jfine + 1][ifine], 2));
						val = val1 + val2 + val3 + val4;

						xx2[ign][k][j][i] = (xx2[ign + 1][kfine][jfine][ifine] * val1 + xx2[ign + 1][kfine][jfine + 1][ifine] * val2
							+ xx2[ign + 1][kfine + 1][jfine][ifine] * val3 + xx2[ign + 1][kfine + 1][jfine + 1][ifine] * val4) / val;
						yy2[ign][k][j][i] = (yy2[ign + 1][kfine][jfine][ifine] * val1 + yy2[ign + 1][kfine][jfine + 1][ifine] * val2
							+ yy2[ign + 1][kfine + 1][jfine][ifine] * val3 + yy2[ign + 1][kfine + 1][jfine + 1][ifine] * val4) / val;
						zz2[ign][k][j][i] = (zz2[ign + 1][kfine][jfine][ifine] * val1 + zz2[ign + 1][kfine][jfine + 1][ifine] * val2
							+ zz2[ign + 1][kfine + 1][jfine][ifine] * val3 + zz2[ign + 1][kfine + 1][jfine + 1][ifine] * val4) / val;



						s2xn[ign][k][j][i] = s2xn[ign + 1][kfine][jfine][ifine] + s2xn[ign + 1][kfine][jfine + 1][ifine]
							+ s2xn[ign + 1][kfine + 1][jfine][ifine] + s2xn[ign + 1][kfine + 1][jfine + 1][ifine];
						s2yn[ign][k][j][i] = s2yn[ign + 1][kfine][jfine][ifine] + s2yn[ign + 1][kfine][jfine + 1][ifine]
							+ s2yn[ign + 1][kfine + 1][jfine][ifine] + s2yn[ign + 1][kfine + 1][jfine + 1][ifine];
						s2zn[ign][k][j][i] = s2zn[ign + 1][kfine][jfine][ifine] + s2zn[ign + 1][kfine][jfine + 1][ifine]
							+ s2zn[ign + 1][kfine + 1][jfine][ifine] + s2zn[ign + 1][kfine + 1][jfine + 1][ifine];
					}
				}
			}



			for (int k = 1; k<nz + 1; k++)
			{

				kfine = 2 * k - 1;
				for (int j = 1; j<ny + 2; j++)
				{
					jfine = 2 * j - 1;
					for (int i = 1; i<nx + 1; i++)
					{

						ifine = 2 * i - 1;

						val1 = sqrt(pow(s3xn[ign + 1][kfine][jfine][ifine], 2) + pow(s3yn[ign + 1][kfine][jfine][ifine], 2) + pow(s3zn[ign + 1][kfine][jfine][ifine], 2));
						val2 = sqrt(pow(s3xn[ign + 1][kfine][jfine][ifine + 1], 2) + pow(s3yn[ign + 1][kfine][jfine][ifine + 1], 2) + pow(s3zn[ign + 1][kfine][jfine][ifine + 1], 2));
						val3 = sqrt(pow(s3xn[ign + 1][kfine + 1][jfine][ifine], 2) + pow(s3yn[ign + 1][kfine + 1][jfine][ifine], 2) + pow(s3zn[ign + 1][kfine + 1][jfine][ifine], 2));
						val4 = sqrt(pow(s3xn[ign + 1][kfine + 1][jfine][ifine + 1], 2) + pow(s3yn[ign + 1][kfine + 1][jfine][ifine + 1], 2) + pow(s3zn[ign + 1][kfine + 1][jfine][ifine + 1], 2));
						val = val1 + val2 + val3 + val4;

						xx3[ign][k][j][i] = (xx3[ign + 1][kfine][jfine][ifine] * val1 + xx3[ign + 1][kfine][jfine][ifine + 1] * val2
							+ xx3[ign + 1][kfine + 1][jfine][ifine] * val3 + xx3[ign + 1][kfine + 1][jfine][ifine + 1] * val4) / val;

						yy3[ign][k][j][i] = (yy3[ign + 1][kfine][jfine][ifine] * val1 + yy3[ign + 1][kfine][jfine][ifine + 1] * val2
							+ yy3[ign + 1][kfine + 1][jfine][ifine] * val3 + yy3[ign + 1][kfine + 1][jfine][ifine + 1] * val4) / val;

						zz3[ign][k][j][i] = (zz3[ign + 1][kfine][jfine][ifine] * val1 + zz3[ign + 1][kfine][jfine][ifine + 1] * val2
							+ zz3[ign + 1][kfine + 1][jfine][ifine] * val3 + zz3[ign + 1][kfine + 1][jfine][ifine + 1] * val4) / val;



						s3xn[ign][k][j][i] = s3xn[ign + 1][kfine][jfine][ifine] + s3xn[ign + 1][kfine][jfine][ifine + 1]
							+ s3xn[ign + 1][kfine + 1][jfine][ifine] + s3xn[ign + 1][kfine + 1][jfine][ifine + 1];

						s3yn[ign][k][j][i] = s3yn[ign + 1][kfine][jfine][ifine] + s3yn[ign + 1][kfine][jfine][ifine + 1]
							+ s3yn[ign + 1][kfine + 1][jfine][ifine] + s3yn[ign + 1][kfine + 1][jfine][ifine + 1];

						s3zn[ign][k][j][i] = s3zn[ign + 1][kfine][jfine][ifine] + s3zn[ign + 1][kfine][jfine][ifine + 1]
							+ s3zn[ign + 1][kfine + 1][jfine][ifine] + s3zn[ign + 1][kfine + 1][jfine][ifine + 1];
					}
				}
			}

			for (int k = 1; k<nz + 2; k++)
			{
				kfine = 2 * k - 1;
				for (int j = 1; j<ny + 1; j++)
				{
					jfine = 2 * j - 1;
					for (int i = 1; i<nx + 1; i++)
					{
						ifine = 2 * i - 1;

						val1 = sqrt(pow(s1xn[ign + 1][kfine][jfine][ifine], 2) + pow(s1yn[ign + 1][kfine][jfine][ifine], 2) + pow(s1zn[ign + 1][kfine][jfine][ifine], 2));
						val2 = sqrt(pow(s1xn[ign + 1][kfine][jfine][ifine + 1], 2) + pow(s1yn[ign + 1][kfine][jfine][ifine + 1], 2) + pow(s1zn[ign + 1][kfine][jfine][ifine + 1], 2));
						val3 = sqrt(pow(s1xn[ign + 1][kfine][jfine + 1][ifine], 2) + pow(s1yn[ign + 1][kfine][jfine + 1][ifine], 2) + pow(s1zn[ign + 1][kfine][jfine + 1][ifine], 2));
						val4 = sqrt(pow(s1xn[ign + 1][kfine][jfine + 1][ifine + 1], 2) + pow(s1yn[ign + 1][kfine][jfine + 1][ifine + 1], 2) + pow(s1zn[ign + 1][kfine][jfine + 1][ifine + 1], 2));
						val = val1 + val2 + val3 + val4;

						xx1[ign][k][j][i] = (xx1[ign + 1][kfine][jfine][ifine] * val1 + xx1[ign + 1][kfine][jfine][ifine + 1] * val2
							+ xx1[ign + 1][kfine][jfine + 1][ifine] * val3 + xx1[ign + 1][kfine][jfine + 1][ifine + 1] * val4) / val;

						yy1[ign][k][j][i] = (yy1[ign + 1][kfine][jfine][ifine] * val1 + yy1[ign + 1][kfine][jfine][ifine + 1] * val2
							+ yy1[ign + 1][kfine][jfine + 1][ifine] * val3 + yy1[ign + 1][kfine][jfine + 1][ifine + 1] * val4) / val;

						zz1[ign][k][j][i] = (zz1[ign + 1][kfine][jfine][ifine] * val1 + zz1[ign + 1][kfine][jfine][ifine + 1] * val2
							+ zz1[ign + 1][kfine][jfine + 1][ifine] * val3 + zz1[ign + 1][kfine][jfine + 1][ifine + 1] * val4) / val;



						s1xn[ign][k][j][i] = s1xn[ign + 1][kfine][jfine][ifine] + s1xn[ign + 1][kfine][jfine][ifine + 1]
							+ s1xn[ign + 1][kfine][jfine + 1][ifine] + s1xn[ign + 1][kfine][jfine + 1][ifine + 1];

						s1yn[ign][k][j][i] = s1yn[ign + 1][kfine][jfine][ifine] + s1yn[ign + 1][kfine][jfine][ifine + 1]
							+ s1yn[ign + 1][kfine][jfine + 1][ifine] + s1yn[ign + 1][kfine][jfine + 1][ifine + 1];

						s1zn[ign][k][j][i] = s1zn[ign + 1][kfine][jfine][ifine] + s1zn[ign + 1][kfine][jfine][ifine + 1]
							+ s1zn[ign + 1][kfine][jfine + 1][ifine] + s1zn[ign + 1][kfine][jfine + 1][ifine + 1];
					}
				}
			}
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					s1xn[ign][0][j][i] = s1xn[ign][1][j][i];
					s1yn[ign][0][j][i] = s1yn[ign][1][j][i];
					s1zn[ign][0][j][i] = s1zn[ign][1][j][i];

					s1xn[ign][nz + 2][j][i] = s1xn[ign][nz + 1][j][i];
					s1yn[ign][nz + 2][j][i] = s1yn[ign][nz + 1][j][i];
					s1zn[ign][nz + 2][j][i] = s1zn[ign][nz + 1][j][i];
				}
			for (int k = 1; k<nz + 1; k++)
				for (int j = 1; j<ny + 1; j++)
				{
					s2xn[ign][k][j][0] = s2xn[ign][k][j][1];
					s2yn[ign][k][j][0] = s2yn[ign][k][j][1];
					s2zn[ign][k][j][0] = s2zn[ign][k][j][1];

					s2xn[ign][k][j][nx + 2] = s2xn[ign][k][j][nx + 1];
					s2yn[ign][k][j][nx + 2] = s2yn[ign][k][j][nx + 1];
					s2zn[ign][k][j][nx + 2] = s2zn[ign][k][j][nx + 1];
				}

			//j-direction 边界值设定


			for (int k = 1; k<nz + 1; k++)
				for (int i = 1; i<nx + 1; i++)
				{
					s3xn[ign][k][0][i] = s3xn[ign][k][1][i];
					s3yn[ign][k][0][i] = s3yn[ign][k][1][i];
					s3zn[ign][k][0][i] = s3zn[ign][k][1][i];

					s3xn[ign][k][ny + 2][i] = s3xn[ign][k][ny + 1][i];
					s3yn[ign][k][ny + 2][i] = s3yn[ign][k][ny + 1][i];
					s3zn[ign][k][ny + 2][i] = s3zn[ign][k][ny + 1][i];
				}
			for (int k = 1; k<nz + 1; k++)
			{
				kfine = 2 * k - 1;
				for (int j = 1; j<ny + 1; j++)
				{
					jfine = 2 * j - 1;
					for (int i = 1; i<nx + 1; i++)
					{
						ifine = 2 * i - 1;

						val1 = vvn[ign + 1][kfine][jfine][ifine];
						val2 = vvn[ign + 1][kfine][jfine][ifine + 1];
						val3 = vvn[ign + 1][kfine][jfine + 1][ifine];
						val4 = vvn[ign + 1][kfine][jfine + 1][ifine + 1];
						val5 = vvn[ign + 1][kfine + 1][jfine][ifine];
						val6 = vvn[ign + 1][kfine + 1][jfine][ifine + 1];
						val7 = vvn[ign + 1][kfine + 1][jfine + 1][ifine];
						val8 = vvn[ign + 1][kfine + 1][jfine + 1][ifine + 1];
						val = val1 + val2 + val3 + val4 + val5 + val6 + val7 + val8;

						xx[ign][k][j][i] = (xx[ign + 1][kfine][jfine][ifine] * val1 + xx[ign + 1][kfine][jfine][ifine + 1] * val2
							+ xx[ign + 1][kfine][jfine + 1][ifine] * val3 + xx[ign + 1][kfine][jfine + 1][ifine + 1] * val4
							+ xx[ign + 1][kfine + 1][jfine][ifine] * val5 + xx[ign + 1][kfine + 1][jfine][ifine + 1] * val6
							+ xx[ign + 1][kfine + 1][jfine + 1][ifine] * val7 + xx[ign + 1][kfine + 1][jfine + 1][ifine + 1] * val8) / val;
						yy[ign][k][j][i] = (yy[ign + 1][kfine][jfine][ifine] * val1 + yy[ign + 1][kfine][jfine][ifine + 1] * val2
							+ yy[ign + 1][kfine][jfine + 1][ifine] * val3 + yy[ign + 1][kfine][jfine + 1][ifine + 1] * val4
							+ yy[ign + 1][kfine + 1][jfine][ifine] * val5 + yy[ign + 1][kfine + 1][jfine][ifine + 1] * val6
							+ yy[ign + 1][kfine + 1][jfine + 1][ifine] * val7 + yy[ign + 1][kfine + 1][jfine + 1][ifine + 1] * val8) / val;

						zz[ign][k][j][i] = (zz[ign + 1][kfine][jfine][ifine] * val1 + zz[ign + 1][kfine][jfine][ifine + 1] * val2
							+ zz[ign + 1][kfine][jfine + 1][ifine] * val3 + zz[ign + 1][kfine][jfine + 1][ifine + 1] * val4
							+ zz[ign + 1][kfine + 1][jfine][ifine] * val5 + zz[ign + 1][kfine + 1][jfine][ifine + 1] * val6
							+ zz[ign + 1][kfine + 1][jfine + 1][ifine] * val7 + zz[ign + 1][kfine + 1][jfine + 1][ifine + 1] * val8) / val;

						vvn[ign][k][j][i] = val;
					}
				}
			}
		}

		temp = 2 / double(lbb)*pi;

		cor1 = cos(temp);
		sir1 = sin(temp);

		for (ign = 1; ign<nng + 1; ign++)
		{
			nx = nnx[ign];
			ny = nny[ign];
			nz = nnz[ign];

			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					yy[ign][nz + 1][j][i] = yy[ign][1][j][i] * cor1 - zz[ign][1][j][i] * sir1;
					zz[ign][nz + 1][j][i] = zz[ign][1][j][i] * cor1 + yy[ign][1][j][i] * sir1;

					yy[ign][0][j][i] = yy[ign][nz][j][i] * cor1 + zz[ign][nz][j][i] * sir1;
					zz[ign][0][j][i] = zz[ign][nz][j][i] * cor1 - yy[ign][nz][j][i] * sir1;
				}
		}
	}

	void geon(int nnng)
	{
		nx = nnx[nnng];
		ny = nny[nnng];
		nz = nnz[nnng];

		ib = nib[nnng];
		it = nit[nnng];
		jb = njb[nnng];
		jt = njt[nnng];

		layer = nnng;


		/*for (int k = 1; k<nz + 2; k++)
		for (int j = 1; j<ny + 1; j++)
		for (int i = 1; i<nx + 1; i++)
		{
		xx01[k][j][i] = xx1[nnng][k][j][i];
		yy01[k][j][i] = yy1[nnng][k][j][i];
		zz01[k][j][i] = zz1[nnng][k][j][i];
		}

		for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 1; j++)
		for (int i = 1; i<nx + 2; i++)
		{
		xx02[k][j][i] = xx2[nnng][k][j][i];
		yy02[k][j][i] = yy2[nnng][k][j][i];
		zz02[k][j][i] = zz2[nnng][k][j][i];
		}

		for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 2; j++)
		for (int i = 1; i<nx + 1; i++)
		{
		xx03[k][j][i] = xx3[nnng][k][j][i];
		yy03[k][j][i] = yy3[nnng][k][j][i];
		zz03[k][j][i] = zz3[nnng][k][j][i];
		}
		//*/
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					xx0[k][j][i] = xx[nnng][k][j][i];
				}

		for (int k = 0; k<nz + 2; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					yy0[k][j][i] = yy[nnng][k][j][i];
					zz0[k][j][i] = zz[nnng][k][j][i];
				}



		/*for (int k = 1; k<nz + 2; k++)
		for (int j = 1; j<ny + 1; j++)
		for (int i = 1; i<nx + 1; i++)
		{
		s1x[k][j][i] = s1xn[nnng][k][j][i];
		s1y[k][j][i] = s1yn[nnng][k][j][i];
		s1z[k][j][i] = s1zn[nnng][k][j][i];
		}

		for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 1; j++)
		for (int i = 1; i<nx + 2; i++)
		{
		s2x[k][j][i] = s2xn[nnng][k][j][i];
		s2y[k][j][i] = s2yn[nnng][k][j][i];
		s2z[k][j][i] = s2zn[nnng][k][j][i];
		}

		for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 2; j++)
		for (int i = 1; i<nx + 1; i++)
		{
		s3x[k][j][i] = s3xn[nnng][k][j][i];
		s3y[k][j][i] = s3yn[nnng][k][j][i];
		s3z[k][j][i] = s3zn[nnng][k][j][i];
		}

		for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 1; j++)
		for (int i = 1; i<nx + 1; i++)
		{
		vv[k][j][i] = vvn[nnng][k][j][i];
		}
		//*/
		//i-direction 边界值设定


		/*for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 1; j++)
		{
		s2x[k][j][0] = s2x[k][j][1];
		s2y[k][j][0] = s2y[k][j][1];
		s2z[k][j][0] = s2z[k][j][1];

		s2x[k][j][nx + 2] = s2x[k][j][nx + 1];
		s2y[k][j][nx + 2] = s2y[k][j][nx + 1];
		s2z[k][j][nx + 2] = s2z[k][j][nx + 1];
		}

		//j-direction 边界值设定


		for (int k = 1; k<nz + 1; k++)
		for (int i = 1; i<nx + 1; i++)
		{
		s3x[k][0][i] = s3x[k][1][i];
		s3y[k][0][i] = s3y[k][1][i];
		s3z[k][0][i] = s3z[k][1][i];

		s3x[k][ny + 2][i] = s3x[k][ny + 1][i];
		s3y[k][ny + 2][i] = s3y[k][ny + 1][i];
		s3z[k][ny + 2][i] = s3z[k][ny + 1][i];
		}

		//k-direction 边界值设定


		for (int j = 1; j<ny + 1; j++)
		for (int i = 1; i<nx + 1; i++)
		{
		s1x[0][j][i] = s1x[1][j][i];
		s1y[0][j][i] = s1y[1][j][i];
		s1z[0][j][i] = s1z[1][j][i];

		s1x[nz + 2][j][i] = s1x[nz + 1][j][i];
		s1y[nz + 2][j][i] = s1y[nz + 1][j][i];
		s1z[nz + 2][j][i] = s1z[nz + 1][j][i];
		}
		//*/
	}
	//边界值


	void ye()
	{
		//double sir,cor;

		if (nng>nng0)
		{
			ny = nny[nng];
			nz = nnz[nng];
		}

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
			{
				y1 = yy2[nng][k][j][1];
				z1 = zz2[nng][k][j][1];

				rr = sqrt(y1*y1 + z1*z1);

				sir = z1 / rr;
				cor = y1 / rr;

				for (int i = 1; i<nt + 1; i++)
				{
					betaxn[i][nng][k][j] = vxx / vee;
					betayn[i][nng][k][j] = (vrr*cor - vtt*sir) / vee;
					betazn[i][nng][k][j] = (vrr*sir + vtt*cor) / vee;

					hatn[i][nng][k][j] = ht;
					petn[i][nng][k][j] = pt;
				}
			}
	}

	void yen(int nnng)
	{
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
			{
				betax[k][j] = betaxn[n][nnng][k][j];
				betay[k][j] = betayn[n][nnng][k][j];
				betaz[k][j] = betazn[n][nnng][k][j];

				hat[k][j] = hatn[n][nnng][k][j];
				pet[k][j] = petn[n][nnng][k][j];
			}
	}

	//时间导数项


	void tsd()
	{

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					ts1[k][j][i] = 0;
					ts2[k][j][i] = 0;
					ts3[k][j][i] = 0;
					ts4[k][j][i] = 0;
					ts5[k][j][i] = 0;
				}

		for (int nn = 1; nn<nt + 1; nn++)
			for (int k = 1; k<nz + 1; k++)
				for (int j = 1; j<ny + 1; j++)
					for (int i = 1; i<nx + 1; i++)
					{
						ts1[k][j][i] = ts1[k][j][i] + dm[nn][n] * q11[nn][k][j][i];
						ts2[k][j][i] = ts2[k][j][i] + dm[nn][n] * q12[nn][k][j][i];
						ts3[k][j][i] = ts3[k][j][i] + dm[nn][n] * q13[nn][k][j][i];
						ts4[k][j][i] = ts4[k][j][i] + dm[nn][n] * q14[nn][k][j][i];
						ts5[k][j][i] = ts5[k][j][i] + dm[nn][n] * q15[nn][k][j][i];
					}
	}

	void tsdsa()
	{
		double temp = 0;
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					ts6[k][j][i] = 0;
				}

		for (int nn = 1; nn<nt + 1; nn++)
			for (int k = 1; k<nz + 1; k++)
				for (int j = 1; j<ny + 1; j++)
					for (int i = 1; i<nx + 1; i++)
					{
						ts6[k][j][i] = ts6[k][j][i] + dm[nn][n] * q16[nn][k][j][i];
					}
	}



	//ppp,bc,step,ddd,qqq,qqqv,pred

	void ppp()
	{
		double qq2, cvl, sir, cor;

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					dim = q11[n][k][j][i];

					pvx[k][j][i] = q12[n][k][j][i] / dim;
					pvy[k][j][i] = q13[n][k][j][i] / dim;
					pvz[k][j][i] = q14[n][k][j][i] / dim;


					vx = pvx[k][j][i];
					vy = pvy[k][j][i];
					vz = pvz[k][j][i];



					y1 = yy0[k][j][i];
					z1 = zz0[k][j][i];

					rr = sqrt(y1*y1 + z1*z1);

					sir = z1 / rr;
					cor = y1 / rr;

					vth[k][j][i] = vz*cor - vy*sir;    //周向速度
					vre[k][j][i] = vz*sir + vy*cor;    //径向速度

					qq2 = vx*vx + vy*vy + vz*vz;

					p[k][j][i] = 0.4*(q15[n][k][j][i] - 0.5*dim*qq2);
					//if(p[k][j][i]<0) cout<<"p minus than 0 in ppp "<<k<<' '<<j<<' '<<i<<' '<<q15[n][k][j][i]<<' '<<dim<<' '<<qq2<<endl;
					t[k][j][i] = p[k][j][i] / (dim*rg);

					cvl = cvl0*pow(t[k][j][i] / t0, 1.5)*(t0 + ts) / (t[k][j][i] + ts);

					q16[n][k][j][i] = max(q16[n][k][j][i], pow(10.0, -4)*cvl);
				}
	}

	//边界条件（周向数据传递），实现并行


	void bc()
	{
		double qq2, sxn, syn, szn, ds, a, deltp, rrhoc;
		double u, v, w, uabs, unorm, rinv, c02, dis, cb, cosa, hb, cc02, cvl;
		double ss, sr, sv, sd, s1, ve, v2, dr;
		double t1, cvlt, cvu, tem, tur1, tur2, tur3, fv1;
		vector<double> hr(ny + 1), hv(ny + 1), hd(ny + 1), hp(ny + 1);
		vector<vector<double> >temp1, temp2;

		//***************杩涘彛杈圭晫***************
		for (int i = 0; i<ny + 1; i++)
		{
			hr[i] = hv[i] = hd[i] = hp[i] = 0;
		}

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
			{
				dim = q11[n][k][j][1];

				vx = q12[n][k][j][1] / dim;
				vy = q13[n][k][j][1] / dim;
				vz = q14[n][k][j][1] / dim;

				qq2 = vx*vx + vy*vy + vz*vz;
				en = q15[n][k][j][1];
				pp = 0.4*(en - 0.5*dim*qq2);
				t1 = pp / (dim*rg);

				cvlt = cvl0*pow(t1 / t0, 1.5)*(t0 + ts) / (t1 + ts);
				cvl = cvlt / dim;
				cvu = c2*cvl;

				tur1 = pow(10.0, -4);
				tur2 = pow(10.0, -6);
				tur3 = pow(10.0, -6);


				while (abs((tur1 - tur3) / tur3)>pow(10.0, -6))
				{
					tem = tur1 / cvl;
					fv1 = 1.0 / (1.0 + pow(cv1 / tem, 3));

					tur2 = cvu / fv1;
					tur3 = tur1;
					tur1 = tur2;
				}

				turi[k][j] = tur2;



				ds = sqrt(s2xn[layer][k][j][1] * s2xn[layer][k][j][1] + s2yn[layer][k][j][1] * s2yn[layer][k][j][1] + s2zn[layer][k][j][1] * s2zn[layer][k][j][1]);

				sxn = s2xn[layer][k][j][1] / ds;
				syn = s2yn[layer][k][j][1] / ds;
				szn = s2zn[layer][k][j][1] / ds;

				u = pvx[k][j][1];
				v = pvy[k][j][1];
				w = pvz[k][j][1];

				uabs = sqrt(u*u + v*v + w*w);
				unorm = u*sxn + v*syn + w*szn;
				a = sqrt(1.4*p[k][j][1] / q11[n][k][j][1]);

				if (uabs<pow(10.0, -20))
					cosa = 1.0;
				else
					cosa = -unorm / uabs;

				rinv = unorm - 5 * a;
				c02 = a*a + 0.2*uabs*uabs;
				dis = (0.4*cosa*cosa + 2.0)*c02 / (0.4*rinv*rinv) - 0.2;

				if (dis<0)
					dis = pow(10.0, -20);

				cb = -rinv*(0.4 / (0.4*cosa*cosa + 2.0))*(1.0 + cosa*sqrt(dis));
				cc02 = min(cb*cb / c02, 1.0);
				hb = hat[k][j] * cc02;

				t[k][j][0] = hb / cp;
				p[k][j][0] = pet[k][j] * pow(cc02, 3.5);
				if (p[k][j][0]<0) cout << "p minus than 0 in bc pre  " << pet[k][j] << ' ' << cc02 << ' ' << pow(cc02, 3.5) << endl;
				q11[n][k][j][0] = p[k][j][0] / (t[k][j][0] * rg);

				uabs = 2.0*(hat[k][j] - hb);

				q15[n][k][j][0] = 2.5*p[k][j][0] + 0.5*q11[n][k][j][0] * uabs;

				pvx[k][j][0] = sqrt(uabs)*betax[k][j];
				pvy[k][j][0] = sqrt(uabs)*betay[k][j];
				pvz[k][j][0] = sqrt(uabs)*betaz[k][j];

				q12[n][k][j][0] = q11[n][k][j][0] * pvx[k][j][0];
				q13[n][k][j][0] = q11[n][k][j][0] * pvy[k][j][0];
				q14[n][k][j][0] = q11[n][k][j][0] * pvz[k][j][0];
				q16[n][k][j][0] = turi[k][j] * q11[n][k][j][0];
			}

		//***************鍑哄彛杈圭晫***************

		for (int j = 1; j<ny + 1; j++)
		{
			ss = 0;
			sr = 0;
			sv = 0;
			sd = 0;

			for (int k = 1; k<nz + 1; k++)
			{
				s1 = sqrt(s2xn[layer][k][j][nx + 1] * s2xn[layer][k][j][nx + 1] + s2yn[layer][k][j][nx + 1] * s2yn[layer][k][j][nx + 1] + s2zn[layer][k][j][nx + 1] * s2zn[layer][k][j][nx + 1]);
				y1 = yy2[layer][k][j][nx + 1];
				z1 = zz2[layer][k][j][nx + 1];

				rr = sqrt(y1*y1 + z1*z1);
				dim = q11[n][k][j][nx];
				ve = vth[k][j][nx];

				ss = ss + s1;
				sr = sr + rr*s1;
				sd = sd + dim*s1;
				sv = sv + ve*s1;
			}

			hr[j] = sr / ss;
			hv[j] = sv / ss;
			hd[j] = sd / ss;
		}

		hp[1] = pb1;

		for (int j = 2; j<ny + 1; j++)
		{
			dim = 0.5*(hd[j - 1] + hd[j]);
			v2 = 0.5*(hv[j - 1] + hv[j]);
			dr = hr[j] - hr[j - 1];
			rr = 0.5*(hr[j] + hr[j - 1]);
			hp[j] = hp[j - 1] + dim*v2*v2 / rr*dr;
		}

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
			{
				peb[k][j] = hp[j];
			}

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
			{
				ds = sqrt(s2xn[layer][k][j][nx + 1] * s2xn[layer][k][j][nx + 1] + s2yn[layer][k][j][nx + 1] * s2yn[layer][k][j][nx + 1] + s2zn[layer][k][j][nx + 1] * s2zn[layer][k][j][nx + 1]);

				sxn = -s2xn[layer][k][j][nx + 1] / ds;
				syn = -s2yn[layer][k][j][nx + 1] / ds;
				szn = -s2zn[layer][k][j][nx + 1] / ds;

				a = sqrt(1.4*p[k][j][nx] / q11[n][k][j][nx]);
				rrhoc = 1.0 / (q11[n][k][j][nx] * a);
				deltp = p[k][j][nx] - peb[k][j];



				p[k][j][nx + 1] = peb[k][j];
				q11[n][k][j][nx + 1] = q11[n][k][j][nx] - deltp / (a*a);

				pvx[k][j][nx + 1] = pvx[k][j][nx] + sxn*deltp*rrhoc;
				pvy[k][j][nx + 1] = pvy[k][j][nx] + syn*deltp*rrhoc;
				pvz[k][j][nx + 1] = pvz[k][j][nx] + szn*deltp*rrhoc;

				q12[n][k][j][nx + 1] = q11[n][k][j][nx + 1] * pvx[k][j][nx + 1];
				q13[n][k][j][nx + 1] = q11[n][k][j][nx + 1] * pvy[k][j][nx + 1];
				q14[n][k][j][nx + 1] = q11[n][k][j][nx + 1] * pvz[k][j][nx + 1];

				qq2 = pvx[k][j][nx + 1] * pvx[k][j][nx + 1] + pvy[k][j][nx + 1] * pvy[k][j][nx + 1] + pvz[k][j][nx + 1] * pvz[k][j][nx + 1];

				q15[n][k][j][nx + 1] = 2.5*p[k][j][nx + 1] + 0.5*q11[n][k][j][nx + 1] * qq2;
				q16[n][k][j][nx + 1] = q16[n][k][j][nx];

				t[k][j][nx + 1] = p[k][j][nx + 1] / (q11[n][k][j][nx + 1] * rg);
			}
		//***************涓婁笅杈圭晫y***************

		for (int k = 1; k<nz + 1; k++)
			for (int i = 1; i<nx + 1; i++)
			{
				//*****鍙舵牴浣嶇疆****

				p[k][0][i] = p[k][1][i];
				t[k][0][i] = t[k][1][i];

				q11[n][k][0][i] = q11[n][k][1][i];

				pvx[k][0][i] = -pvx[k][1][i];
				pvy[k][0][i] = -2.0*rpm*zz3[layer][k][1][i] - pvy[k][1][i];
				pvz[k][0][i] = 2.0*rpm*yy3[layer][k][1][i] - pvz[k][1][i];

				q12[n][k][0][i] = q11[n][k][0][i] * pvx[k][0][i];
				q13[n][k][0][i] = q11[n][k][0][i] * pvy[k][0][i];
				q14[n][k][0][i] = q11[n][k][0][i] * pvz[k][0][i];

				qq2 = pvx[k][0][i] * pvx[k][0][i] + pvy[k][0][i] * pvy[k][0][i] + pvz[k][0][i] * pvz[k][0][i];

				q15[n][k][0][i] = 2.5*p[k][0][i] + 0.5*q11[n][k][0][i] * qq2;
				q16[n][k][0][i] = -q16[n][k][1][i];

				//*****鍙堕《浣嶇疆****

				p[k][ny + 1][i] = p[k][ny][i];
				t[k][ny + 1][i] = t[k][ny][i];

				q11[n][k][ny + 1][i] = q11[n][k][ny][i];

				pvx[k][ny + 1][i] = -pvx[k][ny][i];
				pvy[k][ny + 1][i] = -pvy[k][ny][i];
				pvz[k][ny + 1][i] = -pvz[k][ny][i];

				q12[n][k][ny + 1][i] = q11[n][k][ny + 1][i] * pvx[k][ny + 1][i];
				q13[n][k][ny + 1][i] = q11[n][k][ny + 1][i] * pvy[k][ny + 1][i];
				q14[n][k][ny + 1][i] = q11[n][k][ny + 1][i] * pvz[k][ny + 1][i];

				qq2 = pvx[k][ny + 1][i] * pvx[k][ny + 1][i] + pvy[k][ny + 1][i] * pvy[k][ny + 1][i] + pvz[k][ny + 1][i] * pvz[k][ny + 1][i];

				q15[n][k][ny + 1][i] = 2.5*p[k][ny + 1][i] + 0.5*q11[n][k][ny + 1][i] * qq2;
				q16[n][k][ny + 1][i] = -q16[n][k][ny][i];
			}
		//***************鍓嶅悗鍥哄杈圭晫z鐨勮绠?**************

		for (int j = jb; j<jt + 1; j++)
			for (int i = ib; i<it + 1; i++)
			{
				q11[n][0][j][i] = q11[n][1][j][i];

				pvx[0][j][i] = -pvx[1][j][i];
				pvy[0][j][i] = -2.0*rpm*zz1[layer][1][j][i] - pvy[1][j][i];
				pvz[0][j][i] = 2.0*rpm*yy1[layer][1][j][i] - pvz[1][j][i];

				p[0][j][i] = p[1][j][i];

				q16[n][0][j][i] = -q16[n][1][j][i];

				q11[n][nz + 1][j][i] = q11[n][nz][j][i];

				pvx[nz + 1][j][i] = -pvx[nz][j][i];
				pvy[nz + 1][j][i] = -2.0*rpm*zz1[layer][nz + 1][j][i] - pvy[nz][j][i];
				pvz[nz + 1][j][i] = 2.0*rpm*yy1[layer][nz + 1][j][i] - pvz[nz][j][i];

				p[nz + 1][j][i] = p[nz][j][i];
				q16[n][nz + 1][j][i] = -q16[n][nz][j][i];
			}
		//***********************淇敼鑷宠繖閲?
		//*********鍛ㄦ湡鎬ц竟鐣屼笂铏氭嫙鐨勭墿鐞嗛噺*****************
		if (myid == 0)
		{
			myidl = numprocs - 1;
			myidr = myid + 1;
		}
		else if (myid == numprocs - 1)
		{
			myidl = myid - 1;
			myidr = 0;
		}
		else
		{
			myidl = myid - 1;
			myidr = myid + 1;
		}

		/*	   for (int i = 1; i < ny+1; i++)
		{
		MPI_Sendrecv(&q11[n][1][i][1], ib-1, MPI_DOUBLE, myidl, 20, &q11[n][nz+1][i][1] ,ib-1, MPI_DOUBLE, myidr, 20, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&pvx[1][i][1], ib-1, MPI_DOUBLE, myidl, 21, &pvx[nz+1][i][1] ,ib-1, MPI_DOUBLE, myidr, 21, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&vth[1][i][1], ib-1, MPI_DOUBLE, myidl, 22, &vth[nz+1][i][1] ,ib-1, MPI_DOUBLE, myidr, 22, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&vre[1][i][1], ib-1, MPI_DOUBLE, myidl, 23, &vre[nz+1][i][1] ,ib-1, MPI_DOUBLE, myidr, 23, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&p[1][i][1], ib-1, MPI_DOUBLE, myidl, 24, &p[nz+1][i][1] ,ib-1, MPI_DOUBLE, myidr, 24, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&q16[n][1][i][1], ib-1, MPI_DOUBLE, myidl, 25, &q16[n][nz+1][i][1] ,ib-1, MPI_DOUBLE, myidr, 25, MPI_COMM_WORLD, &status);
		}
		MPI_Request req[12];
		MPI_Status sta[12];

		/*		for(int i=1;i<ny+1;i++)
		{
		MPI_Isend(&q11[n][1][i][1], ib-1, MPI_DOUBLE, myidl, 20, MPI_COMM_WORLD,&req[0]);
		MPI_Isend(&pvx[1][i][1], ib-1, MPI_DOUBLE, myidl, 21, MPI_COMM_WORLD,&req[1]);
		MPI_Isend(&vth[1][i][1], ib-1, MPI_DOUBLE, myidl, 22, MPI_COMM_WORLD,&req[2]);
		MPI_Isend(&vre[1][i][1], ib-1, MPI_DOUBLE, myidl, 23, MPI_COMM_WORLD,&req[3]);
		MPI_Isend(&p[1][i][1], ib-1, MPI_DOUBLE, myidl, 24, MPI_COMM_WORLD,&req[4]);
		MPI_Isend(&q16[n][1][i][1], ib-1, MPI_DOUBLE, myidl, 25, MPI_COMM_WORLD,&req[5]);

		MPI_Irecv(&q11[n][nz+1][i][1],ib-1, MPI_DOUBLE, myidr, 20, MPI_COMM_WORLD,&req[6]);
		MPI_Irecv(&pvx[nz+1][i][1],ib-1, MPI_DOUBLE, myidr, 21, MPI_COMM_WORLD,&req[7]);
		MPI_Irecv(&vth[nz+1][i][1],ib-1, MPI_DOUBLE, myidr, 22, MPI_COMM_WORLD,&req[8]);
		MPI_Irecv(&vre[nz+1][i][1],ib-1, MPI_DOUBLE, myidr, 23, MPI_COMM_WORLD,&req[9]);
		MPI_Irecv(&p[nz+1][i][1], ib-1, MPI_DOUBLE, myidr, 24, MPI_COMM_WORLD,&req[10]);
		MPI_Irecv(&q16[n][nz+1][i][1],ib-1, MPI_DOUBLE, myidr, 25, MPI_COMM_WORLD,&req[11]);

		MPI_Waitall(12,req,sta);
		}
		for (int i = 1; i < ny+1; i++)

		{
		//					if(myid==0) cout<<"bc 20 begin "<<endl;
		MPI_Sendrecv(&q11[n][1][i][it+1], nx-it, MPI_DOUBLE, myidl, 30, &q11[n][nz+1][i][it+1], nx-it, MPI_DOUBLE, myidr, 30, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&pvx[1][i][it+1], nx-it, MPI_DOUBLE, myidl, 31, &pvx[nz+1][i][it+1], nx-it, MPI_DOUBLE, myidr, 31, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&vth[1][i][it+1], nx-it, MPI_DOUBLE, myidl, 32, &vth[nz+1][i][it+1], nx-it, MPI_DOUBLE, myidr, 32, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&vre[1][i][it+1], nx-it, MPI_DOUBLE, myidl, 33, &vre[nz+1][i][it+1], nx-it, MPI_DOUBLE, myidr, 33, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&p[1][i][it+1], nx-it, MPI_DOUBLE, myidl, 34, &p[nz+1][i][it+1], nx-it, MPI_DOUBLE, myidr, 34, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&q16[n][1][i][it+1], nx-it, MPI_DOUBLE, myidl, 35, &q16[n][nz+1][i][it+1], nx-it, MPI_DOUBLE, myidr, 35, MPI_COMM_WORLD, &status);

		}

		*/		/*			  for(int i=1;i<ny+1;i++)
		{
		MPI_Isend(&q11[n][1][i][it+1], nx-it, MPI_DOUBLE, myidl, 30, MPI_COMM_WORLD,&req[0]);
		MPI_Isend(&pvx[1][i][it+1], nx-it, MPI_DOUBLE, myidl, 31, MPI_COMM_WORLD,&req[1]);
		MPI_Isend(&vth[1][i][it+1], nx-it, MPI_DOUBLE, myidl, 32, MPI_COMM_WORLD,&req[2]);
		MPI_Isend(&vre[1][i][it+1], nx-it, MPI_DOUBLE, myidl, 33, MPI_COMM_WORLD,&req[3]);
		MPI_Isend(&p[1][i][it+1], nx-it, MPI_DOUBLE, myidl, 34, MPI_COMM_WORLD,&req[4]);
		MPI_Isend(&q16[n][1][i][it+1], nx-it, MPI_DOUBLE, myidl, 35, MPI_COMM_WORLD,&req[5]);

		MPI_Irecv(&q11[n][nz+1][i][it+1],nx-it, MPI_DOUBLE, myidr, 30, MPI_COMM_WORLD,&req[6]);
		MPI_Irecv(&pvx[nz+1][i][it+1],nx-it, MPI_DOUBLE, myidr, 31, MPI_COMM_WORLD,&req[7]);
		MPI_Irecv(&vth[nz+1][i][it+1],nx-it, MPI_DOUBLE, myidr, 32, MPI_COMM_WORLD,&req[8]);
		MPI_Irecv(&vre[nz+1][i][it+1],nx-it, MPI_DOUBLE, myidr, 33, MPI_COMM_WORLD,&req[9]);
		MPI_Irecv(&p[nz+1][i][it+1], nx-it, MPI_DOUBLE, myidr, 34, MPI_COMM_WORLD,&req[10]);
		MPI_Irecv(&q16[n][nz+1][i][it+1],nx-it, MPI_DOUBLE, myidr, 35, MPI_COMM_WORLD,&req[11]);

		MPI_Waitall(12,req,sta);
		}*/

		/*
		for(int i=jt+1 ;i<ny+1;i++)
		{
		MPI_Sendrecv(&q11[n][1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 40, &q11[n][nz+1][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 40, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&pvx[1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 41, &pvx[nz+1][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 41, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&vth[1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 42, &vth[nz+1][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 42, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&vre[1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 43, &vre[nz+1][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 43, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&p[1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 44, &p[nz+1][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 44, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&q16[n][1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 45, &q16[n][nz+1][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 45, MPI_COMM_WORLD, &status);}

		/*					  for(int i=jt+1;i<ny+1;i++)
		{
		MPI_Isend(&q11[n][1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 40, MPI_COMM_WORLD,&req[0]);
		MPI_Isend(&pvx[1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 41, MPI_COMM_WORLD,&req[1]);
		MPI_Isend(&vth[1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 42, MPI_COMM_WORLD,&req[2]);
		MPI_Isend(&vre[1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 43, MPI_COMM_WORLD,&req[3]);
		MPI_Isend(&p[1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 44, MPI_COMM_WORLD,&req[4]);
		MPI_Isend(&q16[n][1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 45, MPI_COMM_WORLD,&req[5]);

		MPI_Irecv(&q11[n][nz+1][i][ib],it-ib+1, MPI_DOUBLE, myidr, 40, MPI_COMM_WORLD,&req[6]);
		MPI_Irecv(&pvx[nz+1][i][ib],it-ib+1, MPI_DOUBLE, myidr, 41, MPI_COMM_WORLD,&req[7]);
		MPI_Irecv(&vth[nz+1][i][ib],it-ib+1, MPI_DOUBLE, myidr, 42, MPI_COMM_WORLD,&req[8]);
		MPI_Irecv(&vre[nz+1][i][ib],it-ib+1, MPI_DOUBLE, myidr, 43, MPI_COMM_WORLD,&req[9]);
		MPI_Irecv(&p[nz+1][i][ib], it-ib+1, MPI_DOUBLE, myidr, 44, MPI_COMM_WORLD,&req[10]);
		MPI_Irecv(&q16[n][nz+1][i][ib],it-ib+1, MPI_DOUBLE, myidr, 45, MPI_COMM_WORLD,&req[11]);

		MPI_Waitall(12,req,sta);

		}

		*

		for(int i=1;i<ny+1;i++)
		{
		MPI_Sendrecv(&q11[n][nz][i][1], ib-1, MPI_DOUBLE, myidl, 50, &q11[n][0][i][1] ,ib-1, MPI_DOUBLE, myidr, 50, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&pvx[nz][i][1], ib-1, MPI_DOUBLE, myidl, 51, &pvx[0][i][1] ,ib-1, MPI_DOUBLE, myidr, 51, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&vth[nz][i][1], ib-1, MPI_DOUBLE, myidl, 52, &vth[0][i][1] ,ib-1, MPI_DOUBLE, myidr, 52, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&vre[nz][i][1], ib-1, MPI_DOUBLE, myidl, 53, &vre[0][i][1] ,ib-1, MPI_DOUBLE, myidr, 53, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&p[nz][i][1], ib-1, MPI_DOUBLE, myidl, 54, &p[0][i][1] ,ib-1, MPI_DOUBLE, myidr, 54, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&q16[n][nz][i][1], ib-1, MPI_DOUBLE, myidl, 55, &q16[n][0][i][1] ,ib-1, MPI_DOUBLE, myidr, 55, MPI_COMM_WORLD, &status);

		}

		/*	  for(int i=1;i<ny+1;i++)
		{
		MPI_Isend(&q11[n][nz][i][1], ib-1, MPI_DOUBLE, myidl, 50, MPI_COMM_WORLD,&req[0]);
		MPI_Isend(&pvx[nz][i][1], ib-1, MPI_DOUBLE, myidl, 51, MPI_COMM_WORLD,&req[1]);
		MPI_Isend(&vth[nz][i][1], ib-1, MPI_DOUBLE, myidl, 52, MPI_COMM_WORLD,&req[2]);
		MPI_Isend(&vre[nz][i][1], ib-1, MPI_DOUBLE, myidl, 53, MPI_COMM_WORLD,&req[3]);
		MPI_Isend(&p[nz][i][1], ib-1, MPI_DOUBLE, myidl, 54, MPI_COMM_WORLD,&req[4]);
		MPI_Isend(&q16[n][nz][i][1], ib-1, MPI_DOUBLE, myidl, 55, MPI_COMM_WORLD,&req[5]);

		MPI_Irecv(&q11[n][0][i][1],ib-1, MPI_DOUBLE, myidr, 50, MPI_COMM_WORLD,&req[6]);
		MPI_Irecv(&pvx[0][i][1],ib-1, MPI_DOUBLE, myidr, 51, MPI_COMM_WORLD,&req[7]);
		MPI_Irecv(&vth[0][i][1],ib-1, MPI_DOUBLE, myidr, 52, MPI_COMM_WORLD,&req[8]);
		MPI_Irecv(&vre[0][i][1],ib-1, MPI_DOUBLE, myidr, 53, MPI_COMM_WORLD,&req[9]);
		MPI_Irecv(&p[0][i][1], ib-1, MPI_DOUBLE, myidr, 54, MPI_COMM_WORLD,&req[10]);
		MPI_Irecv(&q16[n][0][i][1],ib-1, MPI_DOUBLE, myidr, 55, MPI_COMM_WORLD,&req[11]);

		MPI_Waitall(12,req,sta);

		}*


		for(int i=1;i<ny+1;i++)
		{
		MPI_Sendrecv(&q11[n][nz][i][it+1], nx-it, MPI_DOUBLE, myidl, 60, &q11[n][0][i][it+1] ,nx-it, MPI_DOUBLE, myidr, 60, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&pvx[nz][i][it+1], nx-it, MPI_DOUBLE, myidl, 61, &pvx[0][i][it+1] ,nx-it, MPI_DOUBLE, myidr, 61, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&vth[nz][i][it+1], nx-it, MPI_DOUBLE, myidl, 62, &vth[0][i][it+1] ,nx-it, MPI_DOUBLE, myidr, 62, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&vre[nz][i][it+1], nx-it, MPI_DOUBLE, myidl, 63, &vre[0][i][it+1] ,nx-it, MPI_DOUBLE, myidr, 63, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&p[nz][i][it+1], nx-it, MPI_DOUBLE, myidl, 64, &p[0][i][it+1] ,nx-it, MPI_DOUBLE, myidr, 64, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&q16[n][nz][i][it+1], nx-it, MPI_DOUBLE, myidl, 65, &q16[n][0][1][it+1] ,(nx-it)*ny, MPI_DOUBLE, myidr, 65, MPI_COMM_WORLD, &status);

		}

		/*	  for(int i=1;i<ny+1;i++)
		{
		MPI_Isend(&q11[n][nz][i][it+1], nx-it, MPI_DOUBLE, myidl, 60, MPI_COMM_WORLD,&req[0]);
		MPI_Isend(&pvx[nz][i][it+1], nx-it, MPI_DOUBLE, myidl, 61, MPI_COMM_WORLD,&req[1]);
		MPI_Isend(&vth[nz][i][it+1], nx-it, MPI_DOUBLE, myidl, 62, MPI_COMM_WORLD,&req[2]);
		MPI_Isend(&vre[nz][i][it+1], nx-it, MPI_DOUBLE, myidl, 63, MPI_COMM_WORLD,&req[3]);
		MPI_Isend(&p[nz][i][it+1], nx-it, MPI_DOUBLE, myidl, 64, MPI_COMM_WORLD,&req[4]);
		MPI_Isend(&q16[n][nz][i][it+1], nx-it, MPI_DOUBLE, myidl, 65, MPI_COMM_WORLD,&req[5]);

		MPI_Irecv(&q11[n][0][i][it+1],nx-it, MPI_DOUBLE, myidr, 60, MPI_COMM_WORLD,&req[6]);
		MPI_Irecv(&pvx[0][i][it+1],nx-it, MPI_DOUBLE, myidr, 61, MPI_COMM_WORLD,&req[7]);
		MPI_Irecv(&vth[0][i][it+1],nx-it, MPI_DOUBLE, myidr, 62, MPI_COMM_WORLD,&req[8]);
		MPI_Irecv(&vre[0][i][it+1],nx-it, MPI_DOUBLE, myidr, 63, MPI_COMM_WORLD,&req[9]);
		MPI_Irecv(&p[0][i][it+1], nx-it, MPI_DOUBLE, myidr, 64, MPI_COMM_WORLD,&req[10]);
		MPI_Irecv(&q16[n][0][i][it+1],nx-it, MPI_DOUBLE, myidr, 65, MPI_COMM_WORLD,&req[11]);

		MPI_Waitall(12,req,sta);

		}*


		for(int i=jt+1;i<ny+1;i++)
		{
		MPI_Sendrecv(&q11[n][nz][i][ib], it-ib+1, MPI_DOUBLE, myidl, 70, &q11[n][0][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 70, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&pvx[nz][i][ib], it-ib+1, MPI_DOUBLE, myidl, 71, &pvx[0][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 71, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&vth[nz][i][ib], it-ib+1, MPI_DOUBLE, myidl, 72, &vth[0][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 72, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&vre[nz][i][ib], it-ib+1, MPI_DOUBLE, myidl, 73, &vre[0][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 73, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&p[nz][i][ib], it-ib+1, MPI_DOUBLE, myidl, 74, &p[0][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 74, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&q16[n][nz][i][ib], it-ib+1, MPI_DOUBLE, myidl, 75, &q16[n][0][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 75, MPI_COMM_WORLD, &status);
		}

		/*	  for(int i=jt+1;i<ny+1;i++)
		{
		MPI_Isend(&q11[n][nz][i][ib], it-ib+1, MPI_DOUBLE, myidl, 70, MPI_COMM_WORLD,&req[0]);
		MPI_Isend(&pvx[nz][i][ib], it-ib+1, MPI_DOUBLE, myidl, 71, MPI_COMM_WORLD,&req[1]);
		MPI_Isend(&vth[nz][i][ib], it-ib+1, MPI_DOUBLE, myidl, 72, MPI_COMM_WORLD,&req[2]);
		MPI_Isend(&vre[nz][i][ib], it-ib+1, MPI_DOUBLE, myidl, 73, MPI_COMM_WORLD,&req[3]);
		MPI_Isend(&p[nz][i][ib], it-ib+1, MPI_DOUBLE, myidl, 74, MPI_COMM_WORLD,&req[4]);
		MPI_Isend(&q16[n][nz][i][ib], it-ib+1, MPI_DOUBLE, myidl, 75, MPI_COMM_WORLD,&req[5]);

		MPI_Irecv(&q11[n][0][i][ib],it-ib+1, MPI_DOUBLE, myidr, 70, MPI_COMM_WORLD,&req[6]);
		MPI_Irecv(&pvx[0][i][ib],it-ib+1, MPI_DOUBLE, myidr, 71, MPI_COMM_WORLD,&req[7]);
		MPI_Irecv(&vth[0][i][ib],it-ib+1, MPI_DOUBLE, myidr, 72, MPI_COMM_WORLD,&req[8]);
		MPI_Irecv(&vre[0][i][ib],it-ib+1, MPI_DOUBLE, myidr, 73, MPI_COMM_WORLD,&req[9]);
		MPI_Irecv(&p[0][i][ib], it-ib+1, MPI_DOUBLE, myidr, 74, MPI_COMM_WORLD,&req[10]);
		MPI_Irecv(&q16[n][0][i][ib],it-ib+1, MPI_DOUBLE, myidr, 75, MPI_COMM_WORLD,&req[11]);

		MPI_Waitall(12,req,sta);

		}*/

		temp1.resize(ib - 1);
		temp2.resize(ib - 1);

		for (int i = 0; i < ib - 1; i++)
		{
			temp1[i].resize(ny);
			temp2[i].resize(ny);

			for (int j = 0; j< ny; j++)
			{
				temp1[i][j] = 0;
				temp2[i][j] = 0;
			}
		}

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = q11[n][1][j + 1][i + 1];
			}



		for (int i = 0; i < ib - 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidl, 20, &temp2[i][0], ny, MPI_DOUBLE, myidr, 20, MPI_COMM_WORLD, &status);


		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				q11[n][nz + 1][j + 1][i + 1] = temp2[i][j];
			}



		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = pvx[1][j + 1][i + 1];
			}

		for (int i = 0; i < ib - 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidl, 21, &temp2[i][0], ny, MPI_DOUBLE, myidr, 21, MPI_COMM_WORLD, &status);

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				pvx[nz + 1][j + 1][i + 1] = temp2[i][j];
			}



		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = vth[1][j + 1][i + 1];
			}

		for (int i = 0; i < ib - 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidl, 22, &temp2[i][0], ny, MPI_DOUBLE, myidr, 22, MPI_COMM_WORLD, &status);

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				vth[nz + 1][j + 1][i + 1] = temp2[i][j];
			}

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = vre[1][j + 1][i + 1];
			}

		for (int i = 0; i < ib - 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidl, 23, &temp2[i][0], ny, MPI_DOUBLE, myidr, 23, MPI_COMM_WORLD, &status);

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				vre[nz + 1][j + 1][i + 1] = temp2[i][j];
			}



		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = p[1][j + 1][i + 1];
			}

		for (int i = 0; i < ib - 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidl, 24, &temp2[i][0], ny, MPI_DOUBLE, myidr, 24, MPI_COMM_WORLD, &status);

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				p[nz + 1][j + 1][i + 1] = temp2[i][j];
			}



		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = q16[n][1][j + 1][i + 1];
			}

		for (int i = 0; i < ib - 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidl, 25, &temp2[i][0], ny, MPI_DOUBLE, myidr, 25, MPI_COMM_WORLD, &status);

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				q16[n][nz + 1][j + 1][i + 1] = temp2[i][j];
			}



		//sendrecv 30

		temp1.resize(nx - it);
		temp2.resize(nx - it);

		for (int i = 0; i < nx - it; i++)
		{
			temp1[i].resize(ny);
			temp2[i].resize(ny);

			for (int j = 0; j<ny; j++)
			{
				temp1[i][j] = 0;
				temp2[i][j] = 0;
			}
		}

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = q11[n][1][j + 1][i + it + 1];
			}

		for (int i = 0; i < nx - it; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidl, 30, &temp2[i][0], ny, MPI_DOUBLE, myidr, 30, MPI_COMM_WORLD, &status);

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				q11[n][nz + 1][j + 1][i + it + 1] = temp2[i][j];
			}

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = pvx[1][j + 1][i + it + 1];
			}

		for (int i = 0; i < nx - it; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidl, 31, &temp2[i][0], ny, MPI_DOUBLE, myidr, 31, MPI_COMM_WORLD, &status);

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				pvx[nz + 1][j + 1][i + it + 1] = temp2[i][j];
			}



		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = vth[1][j + 1][i + it + 1];
			}

		for (int i = 0; i < nx - it; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidl, 32, &temp2[i][0], ny, MPI_DOUBLE, myidr, 32, MPI_COMM_WORLD, &status);

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				vth[nz + 1][j + 1][i + it + 1] = temp2[i][j];
			}



		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = vre[1][j + 1][i + it + 1];
			}

		for (int i = 0; i < nx - it; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidl, 33, &temp2[i][0], ny, MPI_DOUBLE, myidr, 33, MPI_COMM_WORLD, &status);

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				vre[nz + 1][j + 1][i + it + 1] = temp2[i][j];
			}



		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = p[1][j + 1][i + it + 1];
			}

		for (int i = 0; i < nx - it; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidl, 34, &temp2[i][0], ny, MPI_DOUBLE, myidr, 34, MPI_COMM_WORLD, &status);

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				p[nz + 1][j + 1][i + it + 1] = temp2[i][j];
			}



		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = q16[n][1][j + 1][i + it + 1];
			}

		for (int i = 0; i < nx - it; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidl, 35, &temp2[i][0], ny, MPI_DOUBLE, myidr, 35, MPI_COMM_WORLD, &status);

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				q16[n][nz + 1][j + 1][i + it + 1] = temp2[i][j];
			}



		//sendrecv 40

		temp1.resize(it - ib + 1);
		temp2.resize(it - ib + 1);

		for (int i = 0; i < it - ib + 1; i++)
		{
			temp1[i].resize(ny - jt);
			temp2[i].resize(ny - jt);

			for (int j = 0; j<ny - jt; j++)
			{
				temp1[i][j] = 0;
				temp2[i][j] = 0;
			}
		}

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				temp1[i][j] = q11[n][1][j + jt + 1][i + ib];
			}

		for (int i = 0; i < it - ib + 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny - jt, MPI_DOUBLE, myidl, 40, &temp2[i][0], ny - jt, MPI_DOUBLE, myidr, 40, MPI_COMM_WORLD, &status);

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				q11[n][nz + 1][j + jt + 1][i + ib] = temp2[i][j];
			}



		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				temp1[i][j] = pvx[1][j + jt + 1][i + ib];

			}

		for (int i = 0; i < it - ib + 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny - jt, MPI_DOUBLE, myidl, 41, &temp2[i][0], ny - jt, MPI_DOUBLE, myidr, 41, MPI_COMM_WORLD, &status);

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				pvx[nz + 1][j + jt + 1][i + ib] = temp2[i][j];
			}

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				temp1[i][j] = vth[1][j + jt + 1][i + ib];
			}

		for (int i = 0; i < it - ib + 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny - jt, MPI_DOUBLE, myidl, 42, &temp2[i][0], ny - jt, MPI_DOUBLE, myidr, 42, MPI_COMM_WORLD, &status);

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				vth[nz + 1][j + jt + 1][i + ib] = temp2[i][j];
			}

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				temp1[i][j] = vre[1][j + jt + 1][i + ib];
			}

		for (int i = 0; i < it - ib + 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny - jt, MPI_DOUBLE, myidl, 43, &temp2[i][0], ny - jt, MPI_DOUBLE, myidr, 43, MPI_COMM_WORLD, &status);

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				vre[nz + 1][j + jt + 1][i + ib] = temp2[i][j];
			}

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				temp1[i][j] = p[1][j + jt + 1][i + ib];
			}

		for (int i = 0; i < it - ib + 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny - jt, MPI_DOUBLE, myidl, 44, &temp2[i][0], ny - jt, MPI_DOUBLE, myidr, 44, MPI_COMM_WORLD, &status);

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				p[nz + 1][j + jt + 1][i + ib] = temp2[i][j];
			}

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				temp1[i][j] = q16[n][1][j + jt + 1][i + ib];
			}

		for (int i = 0; i < it - ib + 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny - jt, MPI_DOUBLE, myidl, 45, &temp2[i][0], ny - jt, MPI_DOUBLE, myidr, 45, MPI_COMM_WORLD, &status);

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				q16[n][nz + 1][j + jt + 1][i + ib] = temp2[i][j];
			}



		//sendrecv 50

		temp1.resize(ib - 1);
		temp2.resize(ib - 1);
		for (int i = 0; i < ib - 1; i++)
		{
			temp1[i].resize(ny);
			temp2[i].resize(ny);

			for (int j = 0; j<ny; j++)
			{
				temp1[i][j] = 0;
				temp2[i][j] = 0;
			}
		}

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = q11[n][nz][j + 1][i + 1];
			}

		for (int i = 0; i < ib - 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidr, 50, &temp2[i][0], ny, MPI_DOUBLE, myidl, 50, MPI_COMM_WORLD, &status);

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				q11[n][0][j + 1][i + 1] = temp2[i][j];
			}

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = pvx[nz][j + 1][i + 1];
			}

		for (int i = 0; i < ib - 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidr, 51, &temp2[i][0], ny, MPI_DOUBLE, myidl, 51, MPI_COMM_WORLD, &status);

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				pvx[0][j + 1][i + 1] = temp2[i][j];
			}

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = vth[nz][j + 1][i + 1];
			}

		for (int i = 0; i < ib - 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidr, 52, &temp2[i][0], ny, MPI_DOUBLE, myidl, 52, MPI_COMM_WORLD, &status);

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				vth[0][j + 1][i + 1] = temp2[i][j];
			}

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = vre[nz][j + 1][i + 1];
			}

		for (int i = 0; i < ib - 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidr, 53, &temp2[i][0], ny, MPI_DOUBLE, myidl, 53, MPI_COMM_WORLD, &status);

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				vre[0][j + 1][i + 1] = temp2[i][j];
			}

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = p[nz][j + 1][i + 1];
			}

		for (int i = 0; i < ib - 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidr, 54, &temp2[i][0], ny, MPI_DOUBLE, myidl, 54, MPI_COMM_WORLD, &status);

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				p[0][j + 1][i + 1] = temp2[i][j];
			}

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = q16[n][nz][j + 1][i + 1];
			}

		for (int i = 0; i < ib - 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidr, 55, &temp2[i][0], ny, MPI_DOUBLE, myidl, 55, MPI_COMM_WORLD, &status);

		for (int i = 0; i < ib - 1; i++)
			for (int j = 0; j < ny; j++)
			{
				q16[n][0][j + 1][i + 1] = temp2[i][j];
			}



		//sendrecv 60

		temp1.resize(nx - it);
		temp2.resize(nx - it);

		for (int i = 0; i < nx - it; i++)
		{
			temp1[i].resize(ny);
			temp2[i].resize(ny);

			for (int j = 0; j<ny; j++)
			{
				temp1[i][j] = 0;
				temp2[i][j] = 0;
			}
		}

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = q11[n][nz][j + 1][i + it + 1];
			}

		for (int i = 0; i < nx - it; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidr, 60, &temp2[i][0], ny, MPI_DOUBLE, myidl, 60, MPI_COMM_WORLD, &status);

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				q11[n][0][j + 1][i + it + 1] = temp2[i][j];
			}

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = pvx[nz][j + 1][i + it + 1];
			}

		for (int i = 0; i < nx - it; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidr, 61, &temp2[i][0], ny, MPI_DOUBLE, myidl, 61, MPI_COMM_WORLD, &status);

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				pvx[0][j + 1][i + it + 1] = temp2[i][j];
			}

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = vth[nz][j + 1][i + it + 1];
			}

		for (int i = 0; i < nx - it; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidr, 62, &temp2[i][0], ny, MPI_DOUBLE, myidl, 62, MPI_COMM_WORLD, &status);

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				vth[0][j + 1][i + it + 1] = temp2[i][j];
			}


		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = vre[nz][j + 1][i + it + 1];
			}

		for (int i = 0; i < nx - it; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidr, 63, &temp2[i][0], ny, MPI_DOUBLE, myidl, 63, MPI_COMM_WORLD, &status);

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				vre[0][j + 1][i + it + 1] = temp2[i][j];
			}


		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = p[nz][j + 1][i + it + 1];
			}

		for (int i = 0; i < nx - it; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidr, 64, &temp2[i][0], ny, MPI_DOUBLE, myidl, 64, MPI_COMM_WORLD, &status);

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				p[0][j + 1][i + it + 1] = temp2[i][j];
			}

		//cout << myid << "64 done"<<endl;

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				temp1[i][j] = q16[n][nz][j + 1][i + it + 1];
			}

		for (int i = 0; i < nx - it; i++)
			MPI_Sendrecv(&temp1[i][0], ny, MPI_DOUBLE, myidr, 65, &temp2[i][0], ny, MPI_DOUBLE, myidl, 65, MPI_COMM_WORLD, &status);

		for (int i = 0; i < nx - it; i++)
			for (int j = 0; j < ny; j++)
			{
				q16[n][0][j + 1][i + it + 1] = temp2[i][j];
			}

		//sendrecv 70

		temp1.resize(it - ib + 1);
		temp2.resize(it - ib + 1);

		for (int i = 0; i < it - ib + 1; i++)
		{
			temp1[i].resize(ny - jt);
			temp2[i].resize(ny - jt);

			for (int j = 0; j<ny - nt; j++)
			{
				temp1[i][j] = 0;
				temp2[i][j] = 0;
			}

		}

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				temp1[i][j] = q11[n][nz][j + jt + 1][i + ib];
			}

		for (int i = 0; i < it - ib + 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny - jt, MPI_DOUBLE, myidr, 70, &temp2[i][0], ny - jt, MPI_DOUBLE, myidl, 70, MPI_COMM_WORLD, &status);

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				q11[n][0][j + jt + 1][i + ib] = temp2[i][j];
			}

		//cout << myid << "70 done"<<endl;

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				temp1[i][j] = pvx[nz][j + jt + 1][i + ib];
			}

		for (int i = 0; i < it - ib + 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny - jt, MPI_DOUBLE, myidr, 71, &temp2[i][0], ny - jt, MPI_DOUBLE, myidl, 71, MPI_COMM_WORLD, &status);

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				pvx[0][j + jt + 1][i + ib] = temp2[i][j];
			}

		//cout << myid << "71 done"<<endl;

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				temp1[i][j] = vth[nz][j + jt + 1][i + ib];
			}

		for (int i = 0; i < it - ib + 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny - jt, MPI_DOUBLE, myidr, 72, &temp2[i][0], ny - jt, MPI_DOUBLE, myidl, 72, MPI_COMM_WORLD, &status);

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				vth[0][j + jt + 1][i + ib] = temp2[i][j];
			}

		//cout << myid << "72 done"<<endl;

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				temp1[i][j] = vre[nz][j + jt + 1][i + ib];
			}

		for (int i = 0; i < it - ib + 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny - jt, MPI_DOUBLE, myidr, 73, &temp2[i][0], ny - jt, MPI_DOUBLE, myidl, 73, MPI_COMM_WORLD, &status);

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				vre[0][j + jt + 1][i + ib] = temp2[i][j];
			}

		//cout << myid << "73 done"<<endl;

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				temp1[i][j] = p[nz][j + jt + 1][i + ib];
			}

		for (int i = 0; i < it - ib + 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny - jt, MPI_DOUBLE, myidr, 74, &temp2[i][0], ny - jt, MPI_DOUBLE, myidl, 74, MPI_COMM_WORLD, &status);

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				p[0][j + jt + 1][i + ib] = temp2[i][j];
			}

		//cout << myid << "74 done"<<endl;

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				temp1[i][j] = q16[n][nz][j + jt + 1][i + ib];
			}

		for (int i = 0; i < it - ib + 1; i++)
			MPI_Sendrecv(&temp1[i][0], ny - jt, MPI_DOUBLE, myidr, 75, &temp2[i][0], ny - jt, MPI_DOUBLE, myidl, 75, MPI_COMM_WORLD, &status);

		for (int i = 0; i < it - ib + 1; i++)
			for (int j = 0; j < ny - jt; j++)
			{
				q16[n][0][j + jt + 1][i + ib] = temp2[i][j];
			}
		//*********鍛ㄥ悜z **********
		for (int j = 1; j<ny + 1; j++)
			for (int i = 1; i<nx + 1; i++)
			{
				if (!((i >= ib) && (i <= it) && (j >= jb) && (j <= jt)))
				{
					y1 = yy0[0][j][i];
					z1 = zz0[0][j][i];
					rr = sqrt(y1*y1 + z1*z1);
					sir = z1 / rr;
					cor = y1 / rr;
					pvy[0][j][i] = vre[0][j][i] * cor - vth[0][j][i] * sir;
					pvz[0][j][i] = vre[0][j][i] * sir + vth[0][j][i] * cor;
					pvy[0][j][i] = vre[0][j][i] * cor - vth[0][j][i] * sir;
					pvz[0][j][i] = vre[0][j][i] * sir + vth[0][j][i] * cor;
					y1 = yy0[nz + 1][j][i];
					z1 = zz0[nz + 1][j][i];
					rr = sqrt(y1*y1 + z1*z1);
					sir = z1 / rr;
					cor = y1 / rr;
					pvy[nz + 1][j][i] = vre[nz + 1][j][i] * cor - vth[nz + 1][j][i] * sir;
					pvz[nz + 1][j][i] = vre[nz + 1][j][i] * sir + vth[nz + 1][j][i] * cor;
				}
				q12[n][0][j][i] = q11[n][0][j][i] * pvx[0][j][i];
				q13[n][0][j][i] = q11[n][0][j][i] * pvy[0][j][i];
				q14[n][0][j][i] = q11[n][0][j][i] * pvz[0][j][i];

				qq2 = pvx[0][j][i] * pvx[0][j][i] + pvy[0][j][i] * pvy[0][j][i] + pvz[0][j][i] * pvz[0][j][i];

				q15[n][0][j][i] = 2.5*p[0][j][i] + 0.5*q11[n][0][j][i] * qq2;

				t[0][j][i] = p[0][j][i] / (q11[n][0][j][i] * rg);
				q12[n][nz + 1][j][i] = q11[n][nz + 1][j][i] * pvx[nz + 1][j][i];

				q13[n][nz + 1][j][i] = q11[n][nz + 1][j][i] * pvy[nz + 1][j][i];

				q14[n][nz + 1][j][i] = q11[n][nz + 1][j][i] * pvz[nz + 1][j][i];

				qq2 = pvx[nz + 1][j][i] * pvx[nz + 1][j][i] + pvy[nz + 1][j][i] * pvy[nz + 1][j][i] + pvz[nz + 1][j][i] * pvz[nz + 1][j][i];

				q15[n][nz + 1][j][i] = 2.5*p[nz + 1][j][i] + 0.5*q11[n][nz + 1][j][i] * qq2;

				t[nz + 1][j][i] = p[nz + 1][j][i] / (q11[n][nz + 1][j][i] * rg);
			}
	}
	//计算粘性系数和导热系数
	void viscosity(double temp, double q6, double &cv, double &kc)
	{
		double cvl, cvt, fv1, tem;

		cvl = cvl0*pow(temp / t0, 1.5)*(t0 + ts) / (temp + ts);
		if (isnan(cvl))   cout << "cvl in viscosity " << cvl0 << ' ' << temp << ' ' << t0 << ' ' << ts << ' ' << temp / t0 << ' ' << pow(temp / t0, 1.5) << ' ' << temp + ts << endl;
		tem = q6 / cvl;
		fv1 = 1.0 / (1.0 + pow(cv1 / tem, 3));
		cvt = q6*fv1;
		if (isnan(cvt))	 cout << "cvt in viscosity " << q6 << ' ' << fv1 << endl;
		cv = cvl + cvt;
		if (isnan(cv)) cout << "cv in viscosity " << cvl << ' ' << cvt << endl;
		kc = cp*(cvl / prl + cvt / prt);
	}

	//当地时间步长

	void step()
	{
		double tc, td, aa, cv, kc;
		vector<vector<vector<double> > > sri1(nz + 2, vector<vector<double> >(ny + 2, vector<double>(nx + 2)));
		vector<vector<vector<double> > > srj1(nz + 2, vector<vector<double> >(ny + 2, vector<double>(nx + 2)));
		vector<vector<vector<double> > > srk1(nz + 2, vector<vector<double> >(ny + 2, vector<double>(nx + 2)));

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{

					// !******x方向

					vx = 0.5*(pvx[k][j][i] + pvx[k][j][i + 1]);
					vy = 0.5*(pvy[k][j][i] + pvy[k][j][i + 1]);
					vz = 0.5*(pvz[k][j][i] + pvz[k][j][i + 1]);

					wx = vx;
					wy = vy + rpm*zz2[layer][k][j][i + 1];
					wz = vz - rpm*yy2[layer][k][j][i + 1];

					tc = abs(wx*s2xn[layer][k][j][i + 1] + wy*s2yn[layer][k][j][i + 1] + wz*s2zn[layer][k][j][i + 1]);

					//******y方向

					vx = 0.5*(pvx[k][j][i] + pvx[k][j + 1][i]);
					vy = 0.5*(pvy[k][j][i] + pvy[k][j + 1][i]);
					vz = 0.5*(pvz[k][j][i] + pvz[k][j + 1][i]);

					wx = vx;
					wy = vy + rpm*zz3[layer][k][j + 1][i];
					wz = vz - rpm*yy3[layer][k][j + 1][i];

					tc = abs(wx*s3xn[layer][k][j + 1][i] + wy*s3yn[layer][k][j + 1][i] + wz*s3zn[layer][k][j + 1][i]) + tc;

					//******z方向

					vx = 0.5*(pvx[k][j][i] + pvx[k + 1][j][i]);
					vy = 0.5*(pvy[k][j][i] + pvy[k + 1][j][i]);
					vz = 0.5*(pvz[k][j][i] + pvz[k + 1][j][i]);

					wx = vx;
					wy = vy + rpm*zz1[layer][k + 1][j][i];
					wz = vz - rpm*yy1[layer][k + 1][j][i];

					tc = abs(wx*s1xn[layer][k + 1][j][i] + wy*s1yn[layer][k + 1][j][i] + wz*s1zn[layer][k + 1][j][i]) + tc;
					//	if(isnan(td))    cout<<"td in step 1 "<<wx<<' '<<wy<<' '<<wz<<endl;

					aa = sqrt(1.4*p[k][j][i] / q11[n][k][j][i]);
					if (isnan(aa)) cout << "aa in step " << p[k][j][i] << ' ' << q11[n][k][j][i] << endl;
					tc = tc + aa*(sqrt(pow(s1xn[layer][k + 1][j][i], 2) + pow(s1yn[layer][k + 1][j][i], 2) + pow(s1zn[layer][k + 1][j][i], 2)) +
						sqrt(pow(s2xn[layer][k][j][i + 1], 2) + pow(s2yn[layer][k][j][i + 1], 2) + pow(s2zn[layer][k][j][i + 1], 2)) +
						sqrt(pow(s3xn[layer][k][j + 1][i], 2) + pow(s3yn[layer][k][j + 1][i], 2) + pow(s3zn[layer][k][j + 1][i], 2)));
					if (isnan(tc)) cout << "tc in step " << aa << ' ' << s1xn[layer][k + 1][j][i] << ' ' << s1yn[layer][k + 1][j][i] << ' ' << s1zn[layer][k + 1][j][i] << ' ' << s2xn[layer][k][j][i + 1] << ' ' << s2yn[layer][k][j][i + 1] << ' ' << s2zn[layer][k][j][i + 1] << ' ' << s3xn[layer][k][j + 1][i] << ' ' << s3yn[layer][k][j + 1][i] << ' ' << s3zn[layer][k][j + 1][i] << endl;



					td = pow(s1xn[layer][k + 1][j][i], 2) + pow(s1yn[layer][k + 1][j][i], 2) + pow(s1zn[layer][k + 1][j][i], 2) +
						pow(s2xn[layer][k][j][i + 1], 2) + pow(s2yn[layer][k][j][i + 1], 2) + pow(s2zn[layer][k][j][i + 1], 2) +
						pow(s3xn[layer][k][j + 1][i], 2) + pow(s3yn[layer][k][j + 1][i], 2) + pow(s3zn[layer][k][j + 1][i], 2) +
						2.0*(abs(s1xn[layer][k + 1][j][i] * s2xn[layer][k][j][i + 1] + s1yn[layer][k + 1][j][i] * s2yn[layer][k][j][i + 1] + s1zn[layer][k + 1][j][i] * s2zn[layer][k][j][i + 1]) +
							abs(s1xn[layer][k + 1][j][i] * s3xn[layer][k][j + 1][i] + s1yn[layer][k + 1][j][i] * s3yn[layer][k][j + 1][i] + s1zn[layer][k + 1][j][i] * s3zn[layer][k][j + 1][i]) +
							abs(s3xn[layer][k][j + 1][i] * s2xn[layer][k][j][i + 1] + s3yn[layer][k][j + 1][i] * s2yn[layer][k][j][i + 1] + s3zn[layer][k][j + 1][i] * s2zn[layer][k][j][i + 1]));
					//		 if(isnan(td))    cout<<"td in step 1 "<<s1x[k+1][j][i]<<' '<<s1y[k+1][j][i]<<' '<<s1z[k+1][j][i]<<' '<<s2x[k][j][i+1]<<' '<<s2y[k][j][i+1]<<' '<<s2z[k][j][i+1]<<' '<<s3x[k][j + 1][i]<<' '<<s3y[k][j + 1][i]<<' '<<s3z[k][j + 1][i]<<endl;

					viscosity(t[k][j][i], q16[n][k][j][i], cv, kc);

					td = 8.0*cv*td / q11[n][k][j][i] / vvn[layer][k][j][i];
					if (isnan(td))	 cout << "td in step " << cv << ' ' << q11[n][k][j][i] << ' ' << vvn[layer][k][j][i] << endl;
					time[k][j][i] = tc + td;
					//						if(isnan(time[k][j][i])) cout<<"time in step "<<tc<<' '<<td<<' '<<i<<' '<<j<<' '<<k<<' '<<nx<<' '<<ny<<' '<<nz<<endl;

					sri[k][j][i] = tc + td;
					//						if(isnan(sri[k][j][i])) cout<<"sri in step main "<<tc<<' '<<td<<' '<<i<<' '<<j<<' '<<k<<' '<<nx<<' '<<ny<<' '<<nz<<endl;
					srj[k][j][i] = tc + td;
					srk[k][j][i] = tc + td;

				}
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					sri1[k][j][i] = sri[k][j][i] * (1 + sqrt(srj[k][j][i] / sri[k][j][i]) + sqrt(srk[k][j][i] / sri[k][j][i]));
					srj1[k][j][i] = srj[k][j][i] * (1 + sqrt(sri[k][j][i] / srj[k][j][i]) + sqrt(srk[k][j][i] / srj[k][j][i]));
					srk1[k][j][i] = srk[k][j][i] * (1 + sqrt(sri[k][j][i] / srk[k][j][i]) + sqrt(srj[k][j][i] / srk[k][j][i]));
				}
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
			{
				sri1[k][j][0] = sri1[k][j][1];
				srj1[k][j][0] = srj1[k][j][1];
				srk1[k][j][0] = srk1[k][j][1];



				sri1[k][j][nx + 1] = sri1[k][j][nx];
				srj1[k][j][nx + 1] = srj1[k][j][nx];
				srk1[k][j][nx + 1] = srk1[k][j][nx];
			}
		for (int k = 1; k<nz + 1; k++)
			for (int i = 1; i<nx + 1; i++)
			{
				sri1[k][0][i] = sri1[k][1][i];
				srj1[k][0][i] = srj1[k][1][i];
				srk1[k][0][i] = srk1[k][1][i];



				sri1[k][ny + 1][i] = sri1[k][ny][i];
				srj1[k][ny + 1][i] = srj1[k][ny][i];
				srk1[k][ny + 1][i] = srk1[k][ny][i];

			}
		for (int j = 1; j<ny + 1; j++)
			for (int i = 1; i<nx + 1; i++)
			{
				sri1[0][j][i] = sri1[1][j][i];
				srj1[0][j][i] = srj1[1][j][i];
				srk1[0][j][i] = srk1[1][j][i];



				sri1[nz + 1][j][i] = sri1[nz][j][i];
				srj1[nz + 1][j][i] = srj1[nz][j][i];
				srk1[nz + 1][j][i] = srk1[nz][j][i];
			}
		for (int k = 0; k<nz + 2; k++)
			for (int j = 0; j<ny + 2; j++)
				for (int i = 0; i<nx + 2; i++)
				{
					sri[k][j][i] = sri1[k][j][i];
					srj[k][j][i] = srj1[k][j][i];
					srk[k][j][i] = srk1[k][j][i];
					if (isnan(sri[k][j][i])) cout << "sri in step " << i << ' ' << j << ' ' << k << ' ' << nx << ' ' << ny << ' ' << nz << endl;
					if (isnan(srj[k][j][i])) cout << "srj in step " << i << ' ' << j << ' ' << k << endl;
					if (isnan(srk[k][j][i])) cout << "srk in step " << i << ' ' << j << ' ' << k << endl;
				}
	}

	//人工粘性

	void ddd()
	{
		double em2, em4;
		vector<double> dp;
		double flu1, flu2, flu3, flu4, flu5, ram;

		for (int k = 0; k<nz + 2; k++)
			for (int j = 0; j<ny + 2; j++)
				for (int i = 0; i<nx + 2; i++)
				{
					av1[k][j][i] = 0;
					av2[k][j][i] = 0;
					av3[k][j][i] = 0;
					av4[k][j][i] = 0;
					av5[k][j][i] = 0;

					q15[n][k][j][i] = q15[n][k][j][i] + p[k][j][i];
				}
		//i-direction

		dp.resize(nx + 2);

		dp[0] = 0;
		dp[nx + 1] = 0;

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
			{
				for (int i = 1; i<nx + 1; i++)
				{
					dp[i] = abs((p[k][j][i + 1] - 2.0*p[k][j][i] + p[k][j][i - 1]) / (p[k][j][i + 1] + 2.0*p[k][j][i] + p[k][j][i - 1]));
				}

				for (int i = 0; i<nx + 1; i++)
				{
					ram = 0.50*(sri[k][j][i] + sri[k][j][i + 1]);
					if (isnan(ram)) cout << "ram error in ddd " << sri[k][j][i] << ' ' << sri[k][j][i + 1] << ' ' << k << ' ' << j << ' ' << i << endl;
					em2 = a2*ram*max(dp[i], dp[i + 1]);
					//		if(isnan(em2)) cout<<"em2 error in ddd "<<a2<<' '<<ram<<' '<<dp[i]<<' '<<dp[i+1]<<endl;
					em4 = a4*ram;
					em4 = max(em4 - em2, 0.0);



					flu1 = em2*(q11[n][k][j][i + 1] - q11[n][k][j][i]);
					//		if(isnan(flu1)) cout<<"flu1 error in ddd-1 "<<em2<<' '<<q11[n][k][j][i+1]<<' '<<q11[n][k][j][i]<<endl;
					flu2 = em2*(q12[n][k][j][i + 1] - q12[n][k][j][i]);
					flu3 = em2*(q13[n][k][j][i + 1] - q13[n][k][j][i]);
					flu4 = em2*(q14[n][k][j][i + 1] - q14[n][k][j][i]);
					flu5 = em2*(q15[n][k][j][i + 1] - q15[n][k][j][i]);



					if ((i >= 1) && (i <= nx - 1))
					{
						flu1 = flu1 + em4*(q11[n][k][j][i - 1] - 3.0*q11[n][k][j][i] + 3.0*q11[n][k][j][i + 1] - q11[n][k][j][i + 2]);
						//			if(isnan(flu1)) cout<<"flu1 error in ddd-2 "<<q11[n][k][j][i]<<' '<<q11[n][k][j][i+1]<<' '<<q11[n][k][j][i+1]<<endl;
						flu2 = flu2 + em4*(q12[n][k][j][i - 1] - 3.0*q12[n][k][j][i] + 3.0*q12[n][k][j][i + 1] - q12[n][k][j][i + 2]);
						flu3 = flu3 + em4*(q13[n][k][j][i - 1] - 3.0*q13[n][k][j][i] + 3.0*q13[n][k][j][i + 1] - q13[n][k][j][i + 2]);
						flu4 = flu4 + em4*(q14[n][k][j][i - 1] - 3.0*q14[n][k][j][i] + 3.0*q14[n][k][j][i + 1] - q14[n][k][j][i + 2]);
						flu5 = flu5 + em4*(q15[n][k][j][i - 1] - 3.0*q15[n][k][j][i] + 3.0*q15[n][k][j][i + 1] - q15[n][k][j][i + 2]);
					}

					av1[k][j][i] = av1[k][j][i] + flu1;
					av2[k][j][i] = av2[k][j][i] + flu2;
					av3[k][j][i] = av3[k][j][i] + flu3;
					av4[k][j][i] = av4[k][j][i] + flu4;
					av5[k][j][i] = av5[k][j][i] + flu5;

					//	if(isnan(av1[k][j][i])) cout<<"av1 error in i dir cord"<<flu1;

					av1[k][j][i + 1] = av1[k][j][i + 1] - flu1;
					av2[k][j][i + 1] = av2[k][j][i + 1] - flu2;
					av3[k][j][i + 1] = av3[k][j][i + 1] - flu3;
					av4[k][j][i + 1] = av4[k][j][i + 1] - flu4;
					av5[k][j][i + 1] = av5[k][j][i + 1] - flu5;
					// if(isnan(av1[k][j][i])) cout<<"av1 error in i dir-2 cord"<<flu1;
				}

			}

		//j-direction

		dp.resize(ny + 2);
		dp[0] = 0;
		dp[ny + 1] = 0;

		for (int k = 1; k<nz + 1; k++)
			for (int j = 0; j<ny + 1; j++)
			{
				for (int i = 1; i<nx + 1; i++)
				{
					if (j != 0)
						dp[j] = abs((p[k][j + 1][i] - 2.0*p[k][j][i] + p[k][j - 1][i]) / (p[k][j + 1][i] + 2.0*p[k][j][i] + p[k][j - 1][i]));
					if (j != ny)
						dp[j + 1] = abs((p[k][j + 2][i] - 2.0*p[k][j + 1][i] + p[k][j][i]) / (p[k][j + 2][i] + 2.0*p[k][j + 1][i] + p[k][j][i]));
					ram = 0.50*(srj[k][j][i] + srj[k][j + 1][i]);
					em2 = a2*ram*max(dp[j], dp[j + 1]);
					em4 = a4*ram;
					em4 = max(em4 - em2, 0.0);

					flu1 = em2*(q11[n][k][j + 1][i] - q11[n][k][j][i]);
					flu2 = em2*(q12[n][k][j + 1][i] - q12[n][k][j][i]);
					flu3 = em2*(q13[n][k][j + 1][i] - q13[n][k][j][i]);
					flu4 = em2*(q14[n][k][j + 1][i] - q14[n][k][j][i]);
					flu5 = em2*(q15[n][k][j + 1][i] - q15[n][k][j][i]);

					if ((j >= 1) && (j <= ny - 1))
					{
						flu1 = flu1 + em4*(q11[n][k][j - 1][i] - 3.0*q11[n][k][j][i] + 3.0*q11[n][k][j + 1][i] - q11[n][k][j + 2][i]);
						flu2 = flu2 + em4*(q12[n][k][j - 1][i] - 3.0*q12[n][k][j][i] + 3.0*q12[n][k][j + 1][i] - q12[n][k][j + 2][i]);
						flu3 = flu3 + em4*(q13[n][k][j - 1][i] - 3.0*q13[n][k][j][i] + 3.0*q13[n][k][j + 1][i] - q13[n][k][j + 2][i]);
						flu4 = flu4 + em4*(q14[n][k][j - 1][i] - 3.0*q14[n][k][j][i] + 3.0*q14[n][k][j + 1][i] - q14[n][k][j + 2][i]);
						flu5 = flu5 + em4*(q15[n][k][j - 1][i] - 3.0*q15[n][k][j][i] + 3.0*q15[n][k][j + 1][i] - q15[n][k][j + 2][i]);
					}

					av1[k][j][i] = av1[k][j][i] + flu1;
					av2[k][j][i] = av2[k][j][i] + flu2;
					av3[k][j][i] = av3[k][j][i] + flu3;
					av4[k][j][i] = av4[k][j][i] + flu4;
					av5[k][j][i] = av5[k][j][i] + flu5;

					//		if(isnan(av1[k][j][i])) cout<<"av1 error in j dir cord"<<flu1;

					av1[k][j + 1][i] = av1[k][j + 1][i] - flu1;
					av2[k][j + 1][i] = av2[k][j + 1][i] - flu2;
					av3[k][j + 1][i] = av3[k][j + 1][i] - flu3;
					av4[k][j + 1][i] = av4[k][j + 1][i] - flu4;
					av5[k][j + 1][i] = av5[k][j + 1][i] - flu5;
					//  if(isnan(av1[k][j][i])) cout<<"av1 error in j dir-2 cord"<<flu1;
				}
			}

		//k-direction

		dp.resize(nz + 2);
		dp[0] = 0;
		dp[nz + 1] = 0;

		for (int k = 0; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
			{
				for (int i = 1; i<nx + 1; i++)
				{
					if (k != 0)
						dp[k] = abs((p[k + 1][j][i] - 2.0*p[k][j][i] + p[k - 1][j][i]) / (p[k + 1][j][i] + 2.0*p[k][j][i] + p[k - 1][j][i]));
					if (k != nz)
						dp[k + 1] = abs((p[k + 2][j][i] - 2.0*p[k + 1][j][i] + p[k][j][i]) / (p[k + 2][j][i] + 2.0*p[k + 1][j][i] + p[k][j][i]));
					ram = 0.50*(srk[k][j][i] + srk[k + 1][j][i]);
					em2 = a2*ram*max(dp[k], dp[k + 1]);
					em4 = a4*ram;
					em4 = max(em4 - em2, 0.0);



					flu1 = em2*(q11[n][k + 1][j][i] - q11[n][k][j][i]);
					flu2 = em2*(q12[n][k + 1][j][i] - q12[n][k][j][i]);
					flu3 = em2*(q13[n][k + 1][j][i] - q13[n][k][j][i]);
					flu4 = em2*(q14[n][k + 1][j][i] - q14[n][k][j][i]);
					flu5 = em2*(q15[n][k + 1][j][i] - q15[n][k][j][i]);



					if ((k >= 1) && (k <= nz - 1))
					{
						flu1 = flu1 + em4*(q11[n][k - 1][j][i] - 3.0*q11[n][k][j][i] + 3.0*q11[n][k + 1][j][i] - q11[n][k + 2][j][i]);
						flu2 = flu2 + em4*(q12[n][k - 1][j][i] - 3.0*q12[n][k][j][i] + 3.0*q12[n][k + 1][j][i] - q12[n][k + 2][j][i]);
						flu3 = flu3 + em4*(q13[n][k - 1][j][i] - 3.0*q13[n][k][j][i] + 3.0*q13[n][k + 1][j][i] - q13[n][k + 2][j][i]);
						flu4 = flu4 + em4*(q14[n][k - 1][j][i] - 3.0*q14[n][k][j][i] + 3.0*q14[n][k + 1][j][i] - q14[n][k + 2][j][i]);
						flu5 = flu5 + em4*(q15[n][k - 1][j][i] - 3.0*q15[n][k][j][i] + 3.0*q15[n][k + 1][j][i] - q15[n][k + 2][j][i]);
					}

					av1[k][j][i] = av1[k][j][i] + flu1;
					av2[k][j][i] = av2[k][j][i] + flu2;
					av3[k][j][i] = av3[k][j][i] + flu3;
					av4[k][j][i] = av4[k][j][i] + flu4;
					av5[k][j][i] = av5[k][j][i] + flu5;

					//			if(isnan(av1[k][j][i])) cout<<"av1 error in k dir cord"<<flu1;

					av1[k + 1][j][i] = av1[k + 1][j][i] - flu1;
					av2[k + 1][j][i] = av2[k + 1][j][i] - flu2;
					av3[k + 1][j][i] = av3[k + 1][j][i] - flu3;
					av4[k + 1][j][i] = av4[k + 1][j][i] - flu4;
					av5[k + 1][j][i] = av5[k + 1][j][i] - flu5;
					//  if(isnan(av1[k][j][i])) cout<<"av1 error in k dir-2 cord"<<flu1;
				}

			}

		for (int k = 0; k<nz + 2; k++)
			for (int j = 0; j<ny + 2; j++)
				for (int i = 0; i<nx + 2; i++)
				{
					q15[n][k][j][i] = q15[n][k][j][i] - p[k][j][i];
				}
	}

	//人工粘性SA

	void dddsa()
	{
		double temp = 0;
		double em2, em4;
		vector<double> dp;
		double flu6, ram;

		for (int k = 0; k<nz + 2; k++)
			for (int j = 0; j<ny + 2; j++)
				for (int i = 0; i<nx + 2; i++)
				{
					av6[k][j][i] = 0;
				}

		//i-direction

		dp.resize(nx + 2);
		dp[0] = 0;
		dp[nx + 1] = 0;


		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
			{

				for (int i = 1; i<nx + 1; i++)
				{
					dp[i] = abs((p[k][j][i + 1] - 2.0*p[k][j][i] + p[k][j][i - 1]) / (p[k][j][i + 1] + 2.0*p[k][j][i] + p[k][j][i - 1]));
				}

				for (int i = 0; i<nx + 1; i++)
				{
					ram = 0.50*(sri[k][j][i] + sri[k][j][i + 1]);
					em2 = a2*ram*max(dp[i], dp[i + 1]);
					em4 = a4*ram;
					em4 = max(em4 - em2, 0.0);

					flu6 = em2*(q16[n][k][j][i + 1] - q16[n][k][j][i]);

					if ((i >= 1) && (i <= nx - 1))
					{
						flu6 = flu6 + em4*(q16[n][k][j][i - 1] - 3.0*q16[n][k][j][i] + 3.0*q16[n][k][j][i + 1] - q16[n][k][j][i + 2]);
					}

					av6[k][j][i] = av6[k][j][i] + flu6;
					av6[k][j][i + 1] = av6[k][j][i + 1] - flu6;

					temp = temp + av6[k][j][i];
				}

			}

		//j-direction

		dp.resize(ny + 2);
		dp[0] = 0;
		dp[ny + 1] = 0;

		for (int k = 1; k<nz + 1; k++)
			for (int i = 1; i<nx + 1; i++)
			{
				for (int j = 1; j<ny + 1; j++)
				{
					dp[j] = abs((p[k][j + 1][i] - 2.0*p[k][j][i] + p[k][j - 1][i]) / (p[k][j + 1][i] + 2.0*p[k][j][i] + p[k][j - 1][i]));
				}

				for (int j = 0; j<ny + 1; j++)
				{
					ram = 0.50*(srj[k][j][i] + srj[k][j + 1][i]);
					em2 = a2*ram*max(dp[j], dp[j + 1]);
					em4 = a4*ram;
					em4 = max(em4 - em2, 0.0);

					flu6 = em2*(q16[n][k][j + 1][i] - q16[n][k][j][i]);

					if ((j >= 1) && (j <= ny - 1))
					{
						flu6 = flu6 + em4*(q16[n][k][j - 1][i] - 3.0*q16[n][k][j][i] + 3.0*q16[n][k][j + 1][i] - q16[n][k][j + 2][i]);
					}

					av6[k][j][i] = av6[k][j][i] + flu6;
					av6[k][j + 1][i] = av6[k][j + 1][i] - flu6;

					temp = temp + av6[k][j][i];
				}
			}

		//k-direction
		temp = 0;

		dp.resize(nz + 2);
		dp[0] = 0;
		dp[nz + 1] = 0;

		for (int i = 1; i<nx + 1; i++)
			for (int j = 1; j<ny + 1; j++)
			{
				for (int k = 1; k<nz + 1; k++)
				{
					dp[k] = abs((p[k + 1][j][i] - 2.0*p[k][j][i] + p[k - 1][j][i]) / (p[k + 1][j][i] + 2.0*p[k][j][i] + p[k - 1][j][i]));
				}

				for (int k = 0; k<nz + 1; k++)
				{
					ram = 0.50*(srk[k][j][i] + srk[k + 1][j][i]);
					em2 = a2*ram*max(dp[k], dp[k + 1]);
					em4 = a4*ram;
					em4 = max(em4 - em2, 0.0);

					flu6 = em2*(q16[n][k + 1][j][i] - q16[n][k][j][i]);

					if ((k >= 1) && (k <= nz - 1))
					{
						flu6 = flu6 + em4*(q16[n][k - 1][j][i] - 3.0*q16[n][k][j][i] + 3.0*q16[n][k + 1][j][i] - q16[n][k + 2][j][i]);
					}

					av6[k][j][i] = av6[k][j][i] + flu6;
					av6[k + 1][j][i] = av6[k + 1][j][i] - flu6;

					temp = temp + av6[k][j][i];
				}
			}
	}

	//粗网格人工粘性

	void dddc()
	{
		double em2;
		double flu1, flu2, flu3, flu4, flu5;

		for (int k = 0; k<nz + 2; k++)
			for (int j = 0; j<ny + 2; j++)
				for (int i = 0; i<nx + 2; i++)
				{

					av1[k][j][i] = 0;
					av2[k][j][i] = 0;
					av3[k][j][i] = 0;
					av4[k][j][i] = 0;
					av5[k][j][i] = 0;

					q15[n][k][j][i] = q15[n][k][j][i] + p[k][j][i];
				}

		//i-direction

		em2 = a2 / 1024.0 * sri[1][1][1];

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 0; i<nx + 1; i++)
				{

					flu1 = em2*(q11[n][k][j][i + 1] - q11[n][k][j][i]);
					flu2 = em2*(q12[n][k][j][i + 1] - q12[n][k][j][i]);
					flu3 = em2*(q13[n][k][j][i + 1] - q13[n][k][j][i]);
					flu4 = em2*(q14[n][k][j][i + 1] - q14[n][k][j][i]);
					flu5 = em2*(q15[n][k][j][i + 1] - q15[n][k][j][i]);

					av1[k][j][i] = av1[k][j][i] + flu1;
					av2[k][j][i] = av2[k][j][i] + flu2;
					av3[k][j][i] = av3[k][j][i] + flu3;
					av4[k][j][i] = av4[k][j][i] + flu4;
					av5[k][j][i] = av5[k][j][i] + flu5;



					av1[k][j][i + 1] = av1[k][j][i + 1] - flu1;
					av2[k][j][i + 1] = av2[k][j][i + 1] - flu2;
					av3[k][j][i + 1] = av3[k][j][i + 1] - flu3;
					av4[k][j][i + 1] = av4[k][j][i + 1] - flu4;
					av5[k][j][i + 1] = av5[k][j][i + 1] - flu5;

				}

		//j-direction

		em2 = a2 / 1024.0 * srj[1][1][1];

		for (int k = 1; k<nz + 1; k++)
			for (int j = 0; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{

					flu1 = em2*(q11[n][k][j + 1][i] - q11[n][k][j][i]);
					flu2 = em2*(q12[n][k][j + 1][i] - q12[n][k][j][i]);
					flu3 = em2*(q13[n][k][j + 1][i] - q13[n][k][j][i]);
					flu4 = em2*(q14[n][k][j + 1][i] - q14[n][k][j][i]);
					flu5 = em2*(q15[n][k][j + 1][i] - q15[n][k][j][i]);

					av1[k][j][i] = av1[k][j][i] + flu1;
					av2[k][j][i] = av2[k][j][i] + flu2;
					av3[k][j][i] = av3[k][j][i] + flu3;
					av4[k][j][i] = av4[k][j][i] + flu4;
					av5[k][j][i] = av5[k][j][i] + flu5;



					av1[k][j + 1][i] = av1[k][j + 1][i] - flu1;
					av2[k][j + 1][i] = av2[k][j + 1][i] - flu2;
					av3[k][j + 1][i] = av3[k][j + 1][i] - flu3;
					av4[k][j + 1][i] = av4[k][j + 1][i] - flu4;
					av5[k][j + 1][i] = av5[k][j + 1][i] - flu5;
				}

		//k-direction

		em2 = a2 / 1024.0 * srk[1][1][1];

		for (int k = 0; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					flu1 = em2*(q11[n][k + 1][j][i] - q11[n][k][j][i]);
					flu2 = em2*(q12[n][k + 1][j][i] - q12[n][k][j][i]);
					flu3 = em2*(q13[n][k + 1][j][i] - q13[n][k][j][i]);
					flu4 = em2*(q14[n][k + 1][j][i] - q14[n][k][j][i]);
					flu5 = em2*(q15[n][k + 1][j][i] - q15[n][k][j][i]);

					av1[k][j][i] = av1[k][j][i] + flu1;
					av2[k][j][i] = av2[k][j][i] + flu2;
					av3[k][j][i] = av3[k][j][i] + flu3;
					av4[k][j][i] = av4[k][j][i] + flu4;
					av5[k][j][i] = av5[k][j][i] + flu5;



					av1[k + 1][j][i] = av1[k + 1][j][i] - flu1;
					av2[k + 1][j][i] = av2[k + 1][j][i] - flu2;
					av3[k + 1][j][i] = av3[k + 1][j][i] - flu3;
					av4[k + 1][j][i] = av4[k + 1][j][i] - flu4;
					av5[k + 1][j][i] = av5[k + 1][j][i] - flu5;
				}

		for (int k = 0; k<nz + 2; k++)
			for (int j = 0; j<ny + 2; j++)
				for (int i = 0; i<nx + 2; i++)
				{
					q15[n][k][j][i] = q15[n][k][j][i] - p[k][j][i];
				}
	}

	//对流通量
	/*修改部分：循环重叠*/
	void qqq()
	{
		double flu1, flu2, flu3, flu4, flu5, vf, rf;
		for (int k = 0; k<nz + 2; k++)
			for (int j = 0; j<ny + 2; j++)
				for (int i = 0; i<nx + 2; i++)
				{
					qc1[k][j][i] = 0;
					qc2[k][j][i] = 0;
					qc3[k][j][i] = 0;
					qc4[k][j][i] = 0;
					qc5[k][j][i] = 0;
				}

		//******************x方向********************
		int i, j, k;
		for (k = 1; k < nz + 1; k++)
		{
			for (j = 1; j < ny + 1; j++)
			{
				for (i = 1; i < nx + 1; i++)
				{
					vx = 0.50*(pvx[k][j][i] + pvx[k][j][i - 1]);
					vy = 0.50*(pvy[k][j][i] + pvy[k][j][i - 1]);
					vz = 0.50*(pvz[k][j][i] + pvz[k][j][i - 1]);

					vf = vx*s2xn[layer][k][j][i] + vy*s2yn[layer][k][j][i] + vz*s2zn[layer][k][j][i];
					rf = rpm*zz2[layer][k][j][i] * s2yn[layer][k][j][i] - rpm*yy2[layer][k][j][i] * s2zn[layer][k][j][i];
					vf = vf + rf;

					dim = 0.50*(q11[n][k][j][i] + q11[n][k][j][i - 1]);
					pp = 0.50*(p[k][j][i] + p[k][j][i - 1]);
					en = 0.50*(q15[n][k][j][i] + q15[n][k][j][i - 1]);

					flu1 = dim*vf;
					flu2 = flu1*vx + pp*s2xn[layer][k][j][i];
					flu3 = flu1*vy + pp*s2yn[layer][k][j][i];
					flu4 = flu1*vz + pp*s2zn[layer][k][j][i];
					flu5 = (en + pp)*vf - pp*rf;

					qc1[k][j][i] = qc1[k][j][i] + flu1;
					qc2[k][j][i] = qc2[k][j][i] + flu2;
					qc3[k][j][i] = qc3[k][j][i] + flu3;
					qc4[k][j][i] = qc4[k][j][i] + flu4;
					qc5[k][j][i] = qc5[k][j][i] + flu5;

					qc1[k][j][i - 1] = qc1[k][j][i - 1] - flu1;
					qc2[k][j][i - 1] = qc2[k][j][i - 1] - flu2;
					qc3[k][j][i - 1] = qc3[k][j][i - 1] - flu3;
					qc4[k][j][i - 1] = qc4[k][j][i - 1] - flu4;
					qc5[k][j][i - 1] = qc5[k][j][i - 1] - flu5;


					/*y方向*/
					vx = 0.50*(pvx[k][j][i] + pvx[k][j - 1][i]);
					vy = 0.50*(pvy[k][j][i] + pvy[k][j - 1][i]);
					vz = 0.50*(pvz[k][j][i] + pvz[k][j - 1][i]);

					vf = vx*s3xn[layer][k][j][i] + vy*s3yn[layer][k][j][i] + vz*s3zn[layer][k][j][i];
					rf = rpm*zz3[layer][k][j][i] * s3yn[layer][k][j][i] - rpm*yy3[layer][k][j][i] * s3zn[layer][k][j][i];
					vf = vf + rf;

					dim = 0.50*(q11[n][k][j][i] + q11[n][k][j - 1][i]);
					pp = 0.50*(p[k][j][i] + p[k][j - 1][i]);
					en = 0.50*(q15[n][k][j][i] + q15[n][k][j - 1][i]);

					flu1 = dim*vf;
					flu2 = flu1*vx + pp*s3xn[layer][k][j][i];
					flu3 = flu1*vy + pp*s3yn[layer][k][j][i];
					flu4 = flu1*vz + pp*s3zn[layer][k][j][i];
					flu5 = (en + pp)*vf - pp*rf;

					qc1[k][j][i] = qc1[k][j][i] + flu1;
					qc2[k][j][i] = qc2[k][j][i] + flu2;
					qc3[k][j][i] = qc3[k][j][i] + flu3;
					qc4[k][j][i] = qc4[k][j][i] + flu4;
					qc5[k][j][i] = qc5[k][j][i] + flu5;

					qc1[k][j - 1][i] = qc1[k][j - 1][i] - flu1;
					qc2[k][j - 1][i] = qc2[k][j - 1][i] - flu2;
					qc3[k][j - 1][i] = qc3[k][j - 1][i] - flu3;
					qc4[k][j - 1][i] = qc4[k][j - 1][i] - flu4;
					qc5[k][j - 1][i] = qc5[k][j - 1][i] - flu5;

					/*z方向*/
					vx = 0.50*(pvx[k][j][i] + pvx[k - 1][j][i]);
					vy = 0.50*(pvy[k][j][i] + pvy[k - 1][j][i]);
					vz = 0.50*(pvz[k][j][i] + pvz[k - 1][j][i]);

					vf = vx*s1xn[layer][k][j][i] + vy*s1yn[layer][k][j][i] + vz*s1zn[layer][k][j][i];
					rf = rpm*zz1[layer][k][j][i] * s1yn[layer][k][j][i] - rpm*yy1[layer][k][j][i] * s1zn[layer][k][j][i];
					vf = vf + rf;

					dim = 0.50*(q11[n][k][j][i] + q11[n][k - 1][j][i]);
					pp = 0.50*(p[k][j][i] + p[k - 1][j][i]);
					en = 0.50*(q15[n][k][j][i] + q15[n][k - 1][j][i]);



					flu1 = dim*vf;
					flu2 = flu1*vx + pp*s1xn[layer][k][j][i];
					flu3 = flu1*vy + pp*s1yn[layer][k][j][i];
					flu4 = flu1*vz + pp*s1zn[layer][k][j][i];
					flu5 = (en + pp)*vf - pp*rf;


					qc1[k][j][i] = qc1[k][j][i] + flu1;
					qc2[k][j][i] = qc2[k][j][i] + flu2;
					qc3[k][j][i] = qc3[k][j][i] + flu3;
					qc4[k][j][i] = qc4[k][j][i] + flu4;
					qc5[k][j][i] = qc5[k][j][i] + flu5;



					qc1[k - 1][j][i] = qc1[k - 1][j][i] - flu1;
					qc2[k - 1][j][i] = qc2[k - 1][j][i] - flu2;
					qc3[k - 1][j][i] = qc3[k - 1][j][i] - flu3;
					qc4[k - 1][j][i] = qc4[k - 1][j][i] - flu4;
					qc5[k - 1][j][i] = qc5[k - 1][j][i] - flu5;
				}
				/*x边界*/
				vx = 0.50*(pvx[k][j][i] + pvx[k][j][i - 1]);
				vy = 0.50*(pvy[k][j][i] + pvy[k][j][i - 1]);
				vz = 0.50*(pvz[k][j][i] + pvz[k][j][i - 1]);

				vf = vx*s2xn[layer][k][j][i] + vy*s2yn[layer][k][j][i] + vz*s2zn[layer][k][j][i];
				rf = rpm*zz2[layer][k][j][i] * s2yn[layer][k][j][i] - rpm*yy2[layer][k][j][i] * s2zn[layer][k][j][i];
				vf = vf + rf;

				dim = 0.50*(q11[n][k][j][i] + q11[n][k][j][i - 1]);
				pp = 0.50*(p[k][j][i] + p[k][j][i - 1]);
				en = 0.50*(q15[n][k][j][i] + q15[n][k][j][i - 1]);

				flu1 = dim*vf;
				flu2 = flu1*vx + pp*s2xn[layer][k][j][i];
				flu3 = flu1*vy + pp*s2yn[layer][k][j][i];
				flu4 = flu1*vz + pp*s2zn[layer][k][j][i];
				flu5 = (en + pp)*vf - pp*rf;

				qc1[k][j][i] = qc1[k][j][i] + flu1;
				qc2[k][j][i] = qc2[k][j][i] + flu2;
				qc3[k][j][i] = qc3[k][j][i] + flu3;
				qc4[k][j][i] = qc4[k][j][i] + flu4;
				qc5[k][j][i] = qc5[k][j][i] + flu5;

				qc1[k][j][i - 1] = qc1[k][j][i - 1] - flu1;
				qc2[k][j][i - 1] = qc2[k][j][i - 1] - flu2;
				qc3[k][j][i - 1] = qc3[k][j][i - 1] - flu3;
				qc4[k][j][i - 1] = qc4[k][j][i - 1] - flu4;
				qc5[k][j][i - 1] = qc5[k][j][i - 1] - flu5;
			}
		}

		//******************y边界方向********************

		for (k = 1; k < nz + 1; k++)
		{
			for (i = 1; i < nx + 1; i++)
			{
				vx = 0.50*(pvx[k][ny + 1][i] + pvx[k][ny + 1 - 1][i]);
				vy = 0.50*(pvy[k][ny + 1][i] + pvy[k][ny + 1 - 1][i]);
				vz = 0.50*(pvz[k][ny + 1][i] + pvz[k][ny + 1 - 1][i]);

				vf = vx*s3xn[layer][k][ny + 1][i] + vy*s3yn[layer][k][ny + 1][i] + vz*s3zn[layer][k][ny + 1][i];
				rf = rpm*zz3[layer][k][ny + 1][i] * s3yn[layer][k][ny + 1][i] - rpm*yy3[layer][k][ny + 1][i] * s3zn[layer][k][ny + 1][i];
				vf = vf + rf;

				dim = 0.50*(q11[n][k][ny + 1][i] + q11[n][k][ny + 1 - 1][i]);
				pp = 0.50*(p[k][ny + 1][i] + p[k][ny + 1 - 1][i]);
				en = 0.50*(q15[n][k][ny + 1][i] + q15[n][k][ny + 1 - 1][i]);



				flu1 = dim*vf;
				flu2 = flu1*vx + pp*s3xn[layer][k][ny + 1][i];
				flu3 = flu1*vy + pp*s3yn[layer][k][ny + 1][i];
				flu4 = flu1*vz + pp*s3zn[layer][k][ny + 1][i];
				flu5 = (en + pp)*vf - pp*rf;



				qc1[k][ny + 1][i] = qc1[k][ny + 1][i] + flu1;
				qc2[k][ny + 1][i] = qc2[k][ny + 1][i] + flu2;
				qc3[k][ny + 1][i] = qc3[k][ny + 1][i] + flu3;
				qc4[k][ny + 1][i] = qc4[k][ny + 1][i] + flu4;
				qc5[k][ny + 1][i] = qc5[k][ny + 1][i] + flu5;



				qc1[k][ny + 1 - 1][i] = qc1[k][ny + 1 - 1][i] - flu1;
				qc2[k][ny + 1 - 1][i] = qc2[k][ny + 1 - 1][i] - flu2;
				qc3[k][ny + 1 - 1][i] = qc3[k][ny + 1 - 1][i] - flu3;
				qc4[k][ny + 1 - 1][i] = qc4[k][ny + 1 - 1][i] - flu4;
				qc5[k][ny + 1 - 1][i] = qc5[k][ny + 1 - 1][i] - flu5;
			}
		}

		//******************z边界方向********************
		for (int j = 1; j<ny + 1; j++)
			for (int i = 1; i<nx + 1; i++)
			{
				vx = 0.50*(pvx[k][j][i] + pvx[k - 1][j][i]);
				vy = 0.50*(pvy[k][j][i] + pvy[k - 1][j][i]);
				vz = 0.50*(pvz[k][j][i] + pvz[k - 1][j][i]);

				vf = vx*s1xn[layer][k][j][i] + vy*s1yn[layer][k][j][i] + vz*s1zn[layer][k][j][i];
				rf = rpm*zz1[layer][k][j][i] * s1yn[layer][k][j][i] - rpm*yy1[layer][k][j][i] * s1zn[layer][k][j][i];
				vf = vf + rf;

				dim = 0.50*(q11[n][k][j][i] + q11[n][k - 1][j][i]);
				pp = 0.50*(p[k][j][i] + p[k - 1][j][i]);
				en = 0.50*(q15[n][k][j][i] + q15[n][k - 1][j][i]);



				flu1 = dim*vf;
				flu2 = flu1*vx + pp*s1xn[layer][k][j][i];
				flu3 = flu1*vy + pp*s1yn[layer][k][j][i];
				flu4 = flu1*vz + pp*s1zn[layer][k][j][i];
				flu5 = (en + pp)*vf - pp*rf;


				qc1[k][j][i] = qc1[k][j][i] + flu1;
				qc2[k][j][i] = qc2[k][j][i] + flu2;
				qc3[k][j][i] = qc3[k][j][i] + flu3;
				qc4[k][j][i] = qc4[k][j][i] + flu4;
				qc5[k][j][i] = qc5[k][j][i] + flu5;



				qc1[k - 1][j][i] = qc1[k - 1][j][i] - flu1;
				qc2[k - 1][j][i] = qc2[k - 1][j][i] - flu2;
				qc3[k - 1][j][i] = qc3[k - 1][j][i] - flu3;
				qc4[k - 1][j][i] = qc4[k - 1][j][i] - flu4;
				qc5[k - 1][j][i] = qc5[k - 1][j][i] - flu5;
			}
	}

	//对流通量SA

	void qqqsa()
	{
		double temp = 0;
		double flu6, tur, vf, rf;

		for (int k = 0; k<nz + 2; k++)
			for (int j = 0; j<ny + 2; j++)
				for (int i = 0; i<nx + 2; i++)
				{
					qc6[k][j][i] = 0;
				}

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 2; i++)
				{
					vx = 0.50*(pvx[k][j][i] + pvx[k][j][i - 1]);
					vy = 0.50*(pvy[k][j][i] + pvy[k][j][i - 1]);
					vz = 0.50*(pvz[k][j][i] + pvz[k][j][i - 1]);

					vf = vx*s2xn[layer][k][j][i] + vy*s2yn[layer][k][j][i] + vz*s2zn[layer][k][j][i];
					rf = rpm*zz2[layer][k][j][i] * s2yn[layer][k][j][i] - rpm*yy2[layer][k][j][i] * s2zn[layer][k][j][i];
					vf = vf + rf;

					tur = 0.50*(q16[n][k][j][i] + q16[n][k][j][i - 1]);

					flu6 = tur*vf;

					qc6[k][j][i] = qc6[k][j][i] + flu6;
					qc6[k][j][i - 1] = qc6[k][j][i - 1] - flu6;
				}

		//******************y方向*******************
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 2; j++)
				for (int i = 1; i<nx + 1; i++)
				{

					vx = 0.50*(pvx[k][j][i] + pvx[k][j - 1][i]);
					vy = 0.50*(pvy[k][j][i] + pvy[k][j - 1][i]);
					vz = 0.50*(pvz[k][j][i] + pvz[k][j - 1][i]);

					vf = vx*s3xn[layer][k][j][i] + vy*s3yn[layer][k][j][i] + vz*s3zn[layer][k][j][i];
					rf = rpm*zz3[layer][k][j][i] * s3yn[layer][k][j][i] - rpm*yy3[layer][k][j][i] * s3zn[layer][k][j][i];
					vf = vf + rf;

					tur = 0.50*(q16[n][k][j][i] + q16[n][k][j - 1][i]);
					flu6 = tur*vf;

					qc6[k][j][i] = qc6[k][j][i] + flu6;
					qc6[k][j - 1][i] = qc6[k][j - 1][i] - flu6;
				}

		//******************z方向********************
		for (int k = 1; k<nz + 2; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					vx = 0.50*(pvx[k][j][i] + pvx[k - 1][j][i]);
					vy = 0.50*(pvy[k][j][i] + pvy[k - 1][j][i]);
					vz = 0.50*(pvz[k][j][i] + pvz[k - 1][j][i]);

					vf = vx*s1xn[layer][k][j][i] + vy*s1yn[layer][k][j][i] + vz*s1zn[layer][k][j][i];
					rf = rpm*zz1[layer][k][j][i] * s1yn[layer][k][j][i] - rpm*yy1[layer][k][j][i] * s1zn[layer][k][j][i];
					vf = vf + rf;

					tur = 0.50*(q16[n][k][j][i] + q16[n][k - 1][j][i]);

					flu6 = tur*vf;
					qc6[k][j][i] = qc6[k][j][i] + flu6;
					qc6[k - 1][j][i] = qc6[k - 1][j][i] - flu6;
				}
	}

	//  ??????????????????????????
	/*
	* ??????????gradsfaceI(s2x, s3x, s1x, pvx, gradfi, 1); gradfi[1-15,1-nx+1,1-ny+1,1-nz+1]
	* ?????ortran????gradsfaceI(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),s1x(1:nx,1:ny,0:nz+2),
	*   pvx(0:nx+1,0:ny+1,0:nz+1),gradfi(1,1:nx+1,0:ny+1,0:nz+1))
	*/
	void gradsfaceI(double*** si, double*** sj, double*** sk,
		double*** q, double**** Idqd, int m)
	{
		double sx, fg, qav, rvol, sx1, sx2, qav1, qav2;

		for (int k = 0; k<nz + 2; k++)
			for (int j = 0; j<ny + 2; j++)
				for (int i = 1; i<nx + 2; i++)
				{
					Idqd[m][k][j][i] = 0;
					//????
				}
		//左右面通量x
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
			{
				for (int i = 1; i<nx + 1; i++)
				{
					sx = 0.50*(si[k][j][i] + si[k][j][i + 1]);
					fg = q[k][j][i] * sx;
					//????
					Idqd[m][k][j][i] = Idqd[m][k][j][i] - fg;
					Idqd[m][k][j][i + 1] = Idqd[m][k][j][i + 1] + fg;
				}
				sx = si[k][j][1];
				fg = 0.50*(q[k][j][0] + q[k][j][1])*sx;
				Idqd[m][k][j][1] = Idqd[m][k][j][1] + fg;
				sx = si[k][j][nx + 1];
				fg = 0.50*(q[k][j][nx] + q[k][j][nx + 1])*sx;
				Idqd[m][k][j][nx + 1] = Idqd[m][k][j][nx + 1] - fg;
			}
		//上下面通量y
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 2; j++)
			{
				for (int i = 2; i<nx + 1; i++)
				{
					sx1 = 0.50*sj[k][j][i];
					qav1 = 0.50*(q[k][j][i] + q[k][j - 1][i]);
					sx2 = 0.50*sj[k][j][i - 1];
					qav2 = 0.50*(q[k][j][i - 1] + q[k][j - 1][i - 1]);
					fg = sx1*qav1 + sx2*qav2;
					Idqd[m][k][j][i] = Idqd[m][k][j][i] + fg;
					Idqd[m][k][j - 1][i] = Idqd[m][k][j - 1][i] - fg;
				}
				sx = sj[k][j][1];
				qav = 0.50*(q[k][j][1] + q[k][j - 1][1]);
				fg = qav*sx;
				Idqd[m][k][j][1] = Idqd[m][k][j][1] + fg;
				Idqd[m][k][j - 1][1] = Idqd[m][k][j - 1][1] - fg;
				sx = sj[k][j][nx];
				qav = 0.50*(q[k][j][nx] + q[k][j - 1][nx]);
				fg = qav*sx;
				Idqd[m][k][j][nx + 1] = Idqd[m][k][j][nx + 1] + fg;
				Idqd[m][k][j - 1][nx + 1] = Idqd[m][k][j - 1][nx + 1] - fg;
			}
		//????????
		for (int k = 1; k<nz + 2; k++)
			for (int j = 1; j<ny + 1; j++)
			{
				for (int i = 2; i<nx + 1; i++)
				{
					sx1 = 0.50*sk[k][j][i];
					qav1 = 0.50*(q[k][j][i] + q[k - 1][j][i]);
					sx2 = 0.50*sk[k][j][i - 1];
					qav2 = 0.50*(q[k][j][i - 1] + q[k - 1][j][i - 1]);
					fg = sx1*qav1 + sx2*qav2;
					Idqd[m][k][j][i] = Idqd[m][k][j][i] + fg;
					Idqd[m][k - 1][j][i] = Idqd[m][k - 1][j][i] - fg;
				}
				sx = sk[k][j][1];
				qav = 0.50*(q[k][j][1] + q[k - 1][j][1]);
				fg = qav*sx;
				Idqd[m][k][j][1] = Idqd[m][k][j][1] + fg;
				Idqd[m][k - 1][j][1] = Idqd[m][k - 1][j][1] - fg;

				sx = sk[k][j][nx];
				qav = 0.50*(q[k][j][nx] + q[k - 1][j][nx]);
				fg = qav*sx;
				Idqd[m][k][j][nx + 1] = Idqd[m][k][j][nx + 1] + fg;
				Idqd[m][k - 1][j][nx + 1] = Idqd[m][k - 1][j][nx + 1] - fg;
			}

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
			{
				for (int i = 2; i<nx + 1; i++)
				{
					rvol = 2.0 / (vvn[layer][k][j][i] + vvn[layer][k][j][i - 1]);
					Idqd[m][k][j][i] = Idqd[m][k][j][i] * rvol;
				}
				rvol = 1.0 / vvn[layer][k][j][1];
				Idqd[m][k][j][1] = Idqd[m][k][j][1] * rvol;
				rvol = 1.0 / vvn[layer][k][j][nx];
				Idqd[m][k][j][nx + 1] = Idqd[m][k][j][nx + 1] * rvol;
			}
	}
	/*
	* ??????????gradsfaceJ(s2x, s3x, s1x, pvx, gradfj, 1);
	* ?????ortran????gradsfaceJ(s2x(0:nx+2,1:ny,1:nz),s3x(1:nx,0:ny+2,1:nz),
	*   s1x(1:nx,1:ny,0:nz+2),pvx(0:nx+1,0:ny+1,0:nz+1),gradfj(1,0:nx+1,1:ny+1,0:nz+1))
	*/
	void gradsfaceJ(double*** si, double*** sj, double*** sk,
		double*** q, double**** Jdqd, int m)
	{
		double sx, fg, qav, rvol, sx1, sx2, qav1, qav2;

		for (int k = 0; k<nz + 2; k++)
			for (int j = 1; j<ny + 2; j++)
				for (int i = 0; i<nx + 2; i++)
				{
					Jdqd[m][k][j][i] = 0;
				}

		for (int k = 1; k<nz + 1; k++)
			for (int i = 1; i<nx + 1; i++)
			{
				for (int j = 1; j<ny + 1; j++)
				{
					sx = 0.50*(sj[k][j][i] + sj[k][j + 1][i]);
					fg = q[k][j][i] * sx;
					Jdqd[m][k][j][i] = Jdqd[m][k][j][i] - fg;
					Jdqd[m][k][j + 1][i] = Jdqd[m][k][j + 1][i] + fg;
				}
				sx = sj[k][1][i];
				fg = 0.50*(q[k][0][i] + q[k][1][i])*sx;
				Jdqd[m][k][1][i] = Jdqd[m][k][1][i] + fg;
				sx = sj[k][ny + 1][i];
				fg = 0.50*(q[k][ny][i] + q[k][ny + 1][i])*sx;
				Jdqd[m][k][ny + 1][i] = Jdqd[m][k][ny + 1][i] - fg;
			}
		//左右面通量x
		for (int k = 1; k<nz + 1; k++)
			for (int i = 1; i<nx + 2; i++)
			{
				for (int j = 2; j<ny + 1; j++)
				{
					sx1 = 0.50*si[k][j][i];
					qav1 = 0.50*(q[k][j][i] + q[k][j][i - 1]);
					sx2 = 0.50*si[k][j - 1][i];
					qav2 = 0.50*(q[k][j - 1][i] + q[k][j - 1][i - 1]);
					fg = sx1*qav1 + sx2*qav2;
					Jdqd[m][k][j][i] = Jdqd[m][k][j][i] + fg;
					Jdqd[m][k][j][i - 1] = Jdqd[m][k][j][i - 1] - fg;
				}
				sx = si[k][1][i];
				qav = 0.50*(q[k][1][i] + q[k][1][i - 1]);
				fg = qav*sx;
				Jdqd[m][k][1][i] = Jdqd[m][k][1][i] + fg;
				Jdqd[m][k][1][i - 1] = Jdqd[m][k][1][i - 1] - fg;
				sx = si[k][ny][i];
				qav = 0.50*(q[k][ny][i] + q[k][ny][i - 1]);
				fg = qav*sx;
				Jdqd[m][k][ny + 1][i] = Jdqd[m][k][ny + 1][i] + fg;
				Jdqd[m][k][ny + 1][i - 1] = Jdqd[m][k][ny + 1][i - 1] - fg;
			}
		//前后面通量z
		for (int k = 1; k<nz + 2; k++)
			for (int i = 1; i<nx + 1; i++)
			{
				for (int j = 2; j<ny + 1; j++)
				{
					sx1 = 0.50*sk[k][j][i];
					qav1 = 0.50*(q[k][j][i] + q[k - 1][j][i]);
					sx2 = 0.50*sk[k][j - 1][i];
					qav2 = 0.50*(q[k][j - 1][i] + q[k - 1][j - 1][i]);
					fg = sx1*qav1 + sx2*qav2;
					Jdqd[m][k][j][i] = Jdqd[m][k][j][i] + fg;
					Jdqd[m][k - 1][j][i] = Jdqd[m][k - 1][j][i] - fg;
				}
				sx = sk[k][1][i];
				qav = 0.50*(q[k][1][i] + q[k - 1][1][i]);
				fg = qav*sx;
				Jdqd[m][k][1][i] = Jdqd[m][k][1][i] + fg;
				Jdqd[m][k - 1][1][i] = Jdqd[m][k - 1][1][i] - fg;
				sx = sk[k][ny][i];
				qav = 0.50*(q[k][ny][i] + q[k - 1][ny][i]);
				fg = qav*sx;
				Jdqd[m][k][ny + 1][i] = Jdqd[m][k][ny + 1][i] + fg;
				Jdqd[m][k - 1][ny + 1][i] = Jdqd[m][k - 1][ny + 1][i] - fg;
			}

		for (int k = 1; k<nz + 1; k++)
			for (int i = 1; i<nx + 1; i++)
			{
				for (int j = 2; j<ny + 1; j++)
				{
					rvol = 2.0 / (vvn[layer][k][j][i] + vvn[layer][k][j - 1][i]);
					Jdqd[m][k][j][i] = Jdqd[m][k][j][i] * rvol;
				}
				rvol = 1.0 / vvn[layer][k][1][i];
				Jdqd[m][k][1][i] = Jdqd[m][k][1][i] * rvol;
				rvol = 1.0 / vvn[layer][k][ny][i];
				Jdqd[m][k][ny + 1][i] = Jdqd[m][k][ny + 1][i] * rvol;
			}
	}
	void gradsfaceK(double*** si, double*** sj, double*** sk,
		double*** q, double**** Kdqd, int m)
	{
		double sx, fg, qav, rvol, sx1, sx2, qav1, qav2;

		for (int k = 1; k<nz + 2; k++)
			for (int j = 0; j<ny + 2; j++)
				for (int i = 0; i<nx + 2; i++)
				{
					Kdqd[m][k][j][i] = 0;
				}

		for (int j = 1; j<ny + 1; j++)
			for (int i = 1; i<nx + 1; i++)
			{
				for (int k = 1; k<nz + 1; k++)
				{
					sx = 0.50*(sk[k][j][i] + sk[k + 1][j][i]);
					fg = q[k][j][i] * sx;
					Kdqd[m][k][j][i] = Kdqd[m][k][j][i] - fg;
					Kdqd[m][k + 1][j][i] = Kdqd[m][k + 1][j][i] + fg;
				}
				sx = sk[1][j][i];
				fg = 0.50*(q[0][j][i] + q[1][j][i])*sx;
				Kdqd[m][1][j][i] = Kdqd[m][1][j][i] + fg;
				sx = sk[nz + 1][j][i];
				fg = 0.50*(q[nz][j][i] + q[nz + 1][j][i])*sx;
				Kdqd[m][nz + 1][j][i] = Kdqd[m][nz + 1][j][i] - fg;
			}
		//	 if(myid==0)cout<<"gradfaceK 1"<<endl;

		for (int j = 1; j<ny + 1; j++)
			for (int i = 1; i<nx + 2; i++)
			{
				for (int k = 2; k<nz + 1; k++)
				{
					sx1 = 0.50*si[k][j][i];
					qav1 = 0.50*(q[k][j][i] + q[k][j][i - 1]);
					sx2 = 0.50*si[k - 1][j][i];
					qav2 = 0.50*(q[k - 1][j][i] + q[k - 1][j][i - 1]);
					fg = sx1*qav1 + sx2*qav2;
					Kdqd[m][k][j][i] = Kdqd[m][k][j][i] + fg;
					Kdqd[m][k][j][i - 1] = Kdqd[m][k][j][i - 1] - fg;
				}
				sx = si[1][j][i];
				qav = 0.50*(q[1][j][i] + q[1][j][i - 1]);
				fg = qav*sx;
				Kdqd[m][1][j][i] = Kdqd[m][1][j][i] + fg;
				Kdqd[m][1][j][i - 1] = Kdqd[m][1][j][i - 1] - fg;
				sx = si[nz][j][i];
				qav = 0.50*(q[nz][j][i] + q[nz][j][i - 1]);
				fg = qav*sx;
				Kdqd[m][nz + 1][j][i] = Kdqd[m][nz + 1][j][i] + fg;
				Kdqd[m][nz + 1][j][i - 1] = Kdqd[m][nz + 1][j][i - 1] - fg;
			}
		//	 if(myid==0)cout<<"gradfaceK 2"<<endl;
		//上下面通量y
		for (int j = 1; j<ny + 2; j++)
			for (int i = 1; i<nx + 1; i++)
			{
				for (int k = 2; k<nz + 1; k++)
				{
					sx1 = 0.50*sj[k][j][i];
					qav1 = 0.50*(q[k][j][i] + q[k][j - 1][i]);
					sx2 = 0.50*sj[k - 1][j][i];
					qav2 = 0.50*(q[k - 1][j][i] + q[k - 1][j - 1][i]);
					fg = sx1*qav1 + sx2*qav2;
					Kdqd[m][k][j][i] = Kdqd[m][k][j][i] + fg;
					Kdqd[m][k][j - 1][i] = Kdqd[m][k][j - 1][i] - fg;
				}
				sx = sj[1][j][i];
				qav = 0.50*(q[1][j][i] + q[1][j - 1][i]);
				fg = qav*sx;
				Kdqd[m][1][j][i] = Kdqd[m][1][j][i] + fg;
				Kdqd[m][1][j - 1][i] = Kdqd[m][1][j - 1][i] - fg;

				sx = sj[nz][j][i];
				qav = 0.50*(q[nz][j][i] + q[nz][j - 1][i]);
				fg = qav*sx;
				Kdqd[m][nz + 1][j][i] = Kdqd[m][nz + 1][j][i] + fg;
				Kdqd[m][nz + 1][j - 1][i] = Kdqd[m][nz + 1][j - 1][i] - fg;
			}
		//	 if(myid==0)cout<<"gradfaceK 3"<<endl;

		for (int j = 1; j<ny + 1; j++)
			for (int i = 1; i<nx + 1; i++)
			{
				for (int k = 2; k<nz + 1; k++)
				{
					rvol = 2.0 / (vvn[layer][k][j][i] + vvn[layer][k - 1][j][i]);
					Kdqd[m][k][j][i] = Kdqd[m][k][j][i] * rvol;
				}
				rvol = 1.0 / vvn[layer][1][j][i];
				Kdqd[m][1][j][i] = Kdqd[m][1][j][i] * rvol;
				rvol = 1.0 / vvn[layer][nz][j][i];
				Kdqd[m][nz + 1][j][i] = Kdqd[m][nz + 1][j][i] * rvol;
			}
		//	 if(myid==0)cout<<"gradfaceK 4"<<endl;

	}
	//界面导数值m
	void gradsface()
	{
		for (int i = 1; i<16; i++)
			for (int l = 0; l<nz + 2; l++)
				for (int k = 0; k<ny + 2; k++)
					for (int j = 1; j<nx + 2; j++)
					{
						gradfi[i][l][k][j] = 0;
					}

		for (int i = 1; i<16; i++)
			for (int l = 0; l<nz + 2; l++)
				for (int k = 1; k<ny + 2; k++)
					for (int j = 0; j<nx + 2; j++)
					{
						gradfj[i][l][k][j] = 0;
					}

		for (int i = 1; i<16; i++)
			for (int l = 1; l<nz + 2; l++)
				for (int k = 0; k<ny + 2; k++)
					for (int j = 0; j<nx + 2; j++)
					{
						gradfk[i][l][k][j] = 0;
					}

		gradsfaceI(s2xn[layer], s3xn[layer], s1xn[layer], pvx, gradfi, 1);
		gradsfaceI(s2yn[layer], s3yn[layer], s1yn[layer], pvx, gradfi, 2);
		gradsfaceI(s2zn[layer], s3zn[layer], s1zn[layer], pvx, gradfi, 3);

		gradsfaceI(s2xn[layer], s3xn[layer], s1xn[layer], pvy, gradfi, 4);
		gradsfaceI(s2yn[layer], s3yn[layer], s1yn[layer], pvy, gradfi, 5);
		gradsfaceI(s2zn[layer], s3zn[layer], s1zn[layer], pvy, gradfi, 6);

		gradsfaceI(s2xn[layer], s3xn[layer], s1xn[layer], pvz, gradfi, 7);
		gradsfaceI(s2yn[layer], s3yn[layer], s1yn[layer], pvz, gradfi, 8);
		gradsfaceI(s2zn[layer], s3zn[layer], s1zn[layer], pvz, gradfi, 9);

		gradsfaceI(s2xn[layer], s3xn[layer], s1xn[layer], t, gradfi, 10);
		gradsfaceI(s2yn[layer], s3yn[layer], s1yn[layer], t, gradfi, 11);
		gradsfaceI(s2zn[layer], s3zn[layer], s1zn[layer], t, gradfi, 12);


		gradsfaceJ(s2xn[layer], s3xn[layer], s1xn[layer], pvx, gradfj, 1);
		gradsfaceJ(s2yn[layer], s3yn[layer], s1yn[layer], pvx, gradfj, 2);
		gradsfaceJ(s2zn[layer], s3zn[layer], s1zn[layer], pvx, gradfj, 3);

		gradsfaceJ(s2xn[layer], s3xn[layer], s1xn[layer], pvy, gradfj, 4);
		gradsfaceJ(s2yn[layer], s3yn[layer], s1yn[layer], pvy, gradfj, 5);
		gradsfaceJ(s2zn[layer], s3zn[layer], s1zn[layer], pvy, gradfj, 6);

		gradsfaceJ(s2xn[layer], s3xn[layer], s1xn[layer], pvz, gradfj, 7);
		gradsfaceJ(s2yn[layer], s3yn[layer], s1yn[layer], pvz, gradfj, 8);
		gradsfaceJ(s2zn[layer], s3zn[layer], s1zn[layer], pvz, gradfj, 9);

		gradsfaceJ(s2xn[layer], s3xn[layer], s1xn[layer], t, gradfj, 10);
		gradsfaceJ(s2yn[layer], s3yn[layer], s1yn[layer], t, gradfj, 11);
		gradsfaceJ(s2zn[layer], s3zn[layer], s1zn[layer], t, gradfj, 12);


		gradsfaceK(s2xn[layer], s3xn[layer], s1xn[layer], pvx, gradfk, 1);
		gradsfaceK(s2yn[layer], s3yn[layer], s1yn[layer], pvx, gradfk, 2);
		gradsfaceK(s2zn[layer], s3zn[layer], s1zn[layer], pvx, gradfk, 3);

		gradsfaceK(s2xn[layer], s3xn[layer], s1xn[layer], pvy, gradfk, 4);
		gradsfaceK(s2yn[layer], s3yn[layer], s1yn[layer], pvy, gradfk, 5);
		gradsfaceK(s2zn[layer], s3zn[layer], s1zn[layer], pvy, gradfk, 6);

		gradsfaceK(s2xn[layer], s3xn[layer], s1xn[layer], pvz, gradfk, 7);
		gradsfaceK(s2yn[layer], s3yn[layer], s1yn[layer], pvz, gradfk, 8);
		gradsfaceK(s2zn[layer], s3zn[layer], s1zn[layer], pvz, gradfk, 9);

		gradsfaceK(s2xn[layer], s3xn[layer], s1xn[layer], t, gradfk, 10);
		gradsfaceK(s2yn[layer], s3yn[layer], s1yn[layer], t, gradfk, 11);
		gradsfaceK(s2zn[layer], s3zn[layer], s1zn[layer], t, gradfk, 12);
	}
	//粘性通量
	void qqqv()
	{
		double flu2, flu3, flu4, flu5, tauxx, tauyy, tauzz, tauxy, tauxz, tauyz, phix, phiy, phiz;
		double two3, uav, vav, wav, q16av, tav, mav, kav;
		double timetest_begin, timetest_end;
		timetest_begin = MPI_Wtime();
		two3 = 2.0 / 3.0;

		gradsface();

		for (int k = 0; k<nz + 2; k++)
			for (int j = 0; j<ny + 2; j++)
				for (int i = 0; i<nx + 2; i++)
				{
					qv2[k][j][i] = 0;
					qv3[k][j][i] = 0;
					qv4[k][j][i] = 0;
					qv5[k][j][i] = 0;
				}

		// i-direction -----------------------------------------------------------------
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 2; i++)
				{
					uav = 0.50*(pvx[k][j][i - 1] + pvx[k][j][i]);
					vav = 0.50*(pvy[k][j][i - 1] + pvy[k][j][i]);
					wav = 0.50*(pvz[k][j][i - 1] + pvz[k][j][i]);
					//????
					q16av = 0.50*(q16[n][k][j][i - 1] + q16[n][k][j][i]);
					tav = 0.50*(t[k][j][i - 1] + t[k][j][i]);

					viscosity(tav, q16av, mav, kav);

					tauxx = two3*mav*(2.0*gradfi[1][k][j][i] - gradfi[5][k][j][i] - gradfi[9][k][j][i]);
					tauyy = two3*mav*(2.0*gradfi[5][k][j][i] - gradfi[1][k][j][i] - gradfi[9][k][j][i]);
					tauzz = two3*mav*(2.0*gradfi[9][k][j][i] - gradfi[1][k][j][i] - gradfi[5][k][j][i]);

					tauxy = mav*(gradfi[2][k][j][i] + gradfi[4][k][j][i]);
					tauxz = mav*(gradfi[3][k][j][i] + gradfi[7][k][j][i]);
					tauyz = mav*(gradfi[6][k][j][i] + gradfi[8][k][j][i]);

					phix = uav*tauxx + vav*tauxy + wav*tauxz + kav*gradfi[10][k][j][i];
					phiy = uav*tauxy + vav*tauyy + wav*tauyz + kav*gradfi[11][k][j][i];
					phiz = uav*tauxz + vav*tauyz + wav*tauzz + kav*gradfi[12][k][j][i];

					flu2 = s2xn[layer][k][j][i] * tauxx + s2yn[layer][k][j][i] * tauxy + s2zn[layer][k][j][i] * tauxz;
					flu3 = s2xn[layer][k][j][i] * tauxy + s2yn[layer][k][j][i] * tauyy + s2zn[layer][k][j][i] * tauyz;
					flu4 = s2xn[layer][k][j][i] * tauxz + s2yn[layer][k][j][i] * tauyz + s2zn[layer][k][j][i] * tauzz;
					flu5 = s2xn[layer][k][j][i] * phix + s2yn[layer][k][j][i] * phiy + s2zn[layer][k][j][i] * phiz;

					qv2[k][j][i] = qv2[k][j][i] + flu2;
					qv3[k][j][i] = qv3[k][j][i] + flu3;
					qv4[k][j][i] = qv4[k][j][i] + flu4;
					qv5[k][j][i] = qv5[k][j][i] + flu5;

					qv2[k][j][i - 1] = qv2[k][j][i - 1] - flu2;
					qv3[k][j][i - 1] = qv3[k][j][i - 1] - flu3;
					qv4[k][j][i - 1] = qv4[k][j][i - 1] - flu4;
					qv5[k][j][i - 1] = qv5[k][j][i - 1] - flu5;
				}
		// j-direction -----------------------------------------------------------------
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 2; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					uav = 0.50*(pvx[k][j - 1][i] + pvx[k][j][i]);
					vav = 0.50*(pvy[k][j - 1][i] + pvy[k][j][i]);
					wav = 0.50*(pvz[k][j - 1][i] + pvz[k][j][i]);

					q16av = 0.50*(q16[n][k][j - 1][i] + q16[n][k][j][i]);
					tav = 0.50*(t[k][j - 1][i] + t[k][j][i]);

					viscosity(tav, q16av, mav, kav);

					tauxx = two3*mav*(2.0*gradfj[1][k][j][i] - gradfj[5][k][j][i] - gradfj[9][k][j][i]);
					tauyy = two3*mav*(2.0*gradfj[5][k][j][i] - gradfj[1][k][j][i] - gradfj[9][k][j][i]);
					tauzz = two3*mav*(2.0*gradfj[9][k][j][i] - gradfj[1][k][j][i] - gradfj[5][k][j][i]);

					tauxy = mav*(gradfj[2][k][j][i] + gradfj[4][k][j][i]);
					tauxz = mav*(gradfj[3][k][j][i] + gradfj[7][k][j][i]);
					tauyz = mav*(gradfj[6][k][j][i] + gradfj[8][k][j][i]);

					phix = uav*tauxx + vav*tauxy + wav*tauxz + kav*gradfj[10][k][j][i];
					phiy = uav*tauxy + vav*tauyy + wav*tauyz + kav*gradfj[11][k][j][i];
					phiz = uav*tauxz + vav*tauyz + wav*tauzz + kav*gradfj[12][k][j][i];

					flu2 = s3xn[layer][k][j][i] * tauxx + s3yn[layer][k][j][i] * tauxy + s3zn[layer][k][j][i] * tauxz;
					flu3 = s3xn[layer][k][j][i] * tauxy + s3yn[layer][k][j][i] * tauyy + s3zn[layer][k][j][i] * tauyz;
					flu4 = s3xn[layer][k][j][i] * tauxz + s3yn[layer][k][j][i] * tauyz + s3zn[layer][k][j][i] * tauzz;
					flu5 = s3xn[layer][k][j][i] * phix + s3yn[layer][k][j][i] * phiy + s3zn[layer][k][j][i] * phiz;

					qv2[k][j][i] = qv2[k][j][i] + flu2;
					qv3[k][j][i] = qv3[k][j][i] + flu3;
					qv4[k][j][i] = qv4[k][j][i] + flu4;
					qv5[k][j][i] = qv5[k][j][i] + flu5;

					qv2[k][j - 1][i] = qv2[k][j - 1][i] - flu2;
					qv3[k][j - 1][i] = qv3[k][j - 1][i] - flu3;
					qv4[k][j - 1][i] = qv4[k][j - 1][i] - flu4;
					qv5[k][j - 1][i] = qv5[k][j - 1][i] - flu5;
				}
		// k-direction -----------------------------------------------------------------
		for (int k = 1; k<nz + 2; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					uav = 0.50*(pvx[k - 1][j][i] + pvx[k][j][i]);
					vav = 0.50*(pvy[k - 1][j][i] + pvy[k][j][i]);
					wav = 0.50*(pvz[k - 1][j][i] + pvz[k][j][i]);

					q16av = 0.50*(q16[n][k - 1][j][i] + q16[n][k][j][i]);
					tav = 0.50*(t[k - 1][j][i] + t[k][j][i]);

					viscosity(tav, q16av, mav, kav);

					tauxx = two3*mav*(2.0*gradfk[1][k][j][i] - gradfk[5][k][j][i] - gradfk[9][k][j][i]);
					tauyy = two3*mav*(2.0*gradfk[5][k][j][i] - gradfk[1][k][j][i] - gradfk[9][k][j][i]);
					tauzz = two3*mav*(2.0*gradfk[9][k][j][i] - gradfk[1][k][j][i] - gradfk[5][k][j][i]);

					tauxy = mav*(gradfk[2][k][j][i] + gradfk[4][k][j][i]);
					tauxz = mav*(gradfk[3][k][j][i] + gradfk[7][k][j][i]);
					tauyz = mav*(gradfk[6][k][j][i] + gradfk[8][k][j][i]);

					phix = uav*tauxx + vav*tauxy + wav*tauxz + kav*gradfk[10][k][j][i];
					phiy = uav*tauxy + vav*tauyy + wav*tauyz + kav*gradfk[11][k][j][i];
					phiz = uav*tauxz + vav*tauyz + wav*tauzz + kav*gradfk[12][k][j][i];

					flu2 = s1xn[layer][k][j][i] * tauxx + s1yn[layer][k][j][i] * tauxy + s1zn[layer][k][j][i] * tauxz;
					flu3 = s1xn[layer][k][j][i] * tauxy + s1yn[layer][k][j][i] * tauyy + s1zn[layer][k][j][i] * tauyz;
					flu4 = s1xn[layer][k][j][i] * tauxz + s1yn[layer][k][j][i] * tauyz + s1zn[layer][k][j][i] * tauzz;
					flu5 = s1xn[layer][k][j][i] * phix + s1yn[layer][k][j][i] * phiy + s1zn[layer][k][j][i] * phiz;

					qv2[k][j][i] = qv2[k][j][i] + flu2;
					qv3[k][j][i] = qv3[k][j][i] + flu3;
					qv4[k][j][i] = qv4[k][j][i] + flu4;
					qv5[k][j][i] = qv5[k][j][i] + flu5;

					qv2[k - 1][j][i] = qv2[k - 1][j][i] - flu2;
					qv3[k - 1][j][i] = qv3[k - 1][j][i] - flu3;
					qv4[k - 1][j][i] = qv4[k - 1][j][i] - flu4;
					qv5[k - 1][j][i] = qv5[k - 1][j][i] - flu5;
				}

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					// ????
					qv3[k][j][i] = qv3[k][j][i] + rpm*q14[n][k][j][i] * vvn[layer][k][j][i];
					qv4[k][j][i] = qv4[k][j][i] - rpm*q13[n][k][j][i] * vvn[layer][k][j][i];
				}
		timetest_end = MPI_Wtime();
		if (myid == 0)
		{
			cout << endl;
			cout << "************" << endl;
			cout << " timetest_end-timetest_begin " << timetest_end - timetest_begin << endl;
			cout << "************" << endl;
		}
	}
	//计算方程残差
	void rrr()
	{
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					rr1[k][j][i] = 0;
					rr2[k][j][i] = 0;
					rr3[k][j][i] = 0;
					rr4[k][j][i] = 0;
					rr5[k][j][i] = 0;
				}
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					rr1[k][j][i] = qp1[k][j][i] - qc1[k][j][i] + av1[k][j][i] - ts1[k][j][i] * vvn[layer][k][j][i];
					rr2[k][j][i] = qp2[k][j][i] - qc2[k][j][i] + av2[k][j][i] - ts2[k][j][i] * vvn[layer][k][j][i] + qv2[k][j][i];
					rr3[k][j][i] = qp3[k][j][i] - qc3[k][j][i] + av3[k][j][i] - ts3[k][j][i] * vvn[layer][k][j][i] + qv3[k][j][i];
					rr4[k][j][i] = qp4[k][j][i] - qc4[k][j][i] + av4[k][j][i] - ts4[k][j][i] * vvn[layer][k][j][i] + qv4[k][j][i];
					rr5[k][j][i] = qp5[k][j][i] - qc5[k][j][i] + av5[k][j][i] - ts5[k][j][i] * vvn[layer][k][j][i] + qv5[k][j][i];
					if (isnan(rr1[k][j][i])) cout << "rr1 error in rrr cord" << k << ' ' << j << ' ' << i << ' ' << qp1[k][j][i] << ' ' << qc1[k][j][i] << ' ' << av1[k][j][i] << ' ' << ts1[k][j][i] << ' ' << vvn[layer][k][j][i] << endl;
				}
	}
	//TDMA算法，是一维三节点格式形成的代数方程组的直接求解算法，其实质上是一种特殊的消元-回代法。
	void tdma(vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, vector<double>& x, int n)
	{
		vector<double> u(n + 1), l(n + 1), y(n + 1);
		u[1] = b[1];
		y[1] = d[1];
		for (int i = 2; i<n + 1; i++)
		{
			l[i] = a[i] / u[i - 1];
			u[i] = b[i] - l[i] * c[i - 1];
			y[i] = d[i] - l[i] * y[i - 1];
		}
		x[n] = y[n] / u[n];
		for (int i = n - 1; i>0; i--)
		{
			x[i] = (y[i] - c[i] * x[i + 1]) / u[i];
		}
	}
	//隐式残差光顺
	void ave()
	{
		vector<double> ax(nx + 1), bx(nx + 1), ay(ny + 1), by(ny + 1), az(nz + 1), bz(nz + 1), c1, c2, c3, c4, c5;
		for (int i = 1; i<nx + 1; i++)
		{
			ax[i] = -ta;
			bx[i] = 1.0 + 2.0 * ta;
		}
		for (int i = 1; i<ny + 1; i++)
		{
			ay[i] = -ta;
			by[i] = 1.0 + 2.0 * ta;
		}
		for (int i = 1; i<nz + 1; i++)
		{
			az[i] = -ta;
			bz[i] = 1.0 + 2.0 * ta;
		}
		//????
		c1.resize(nx + 1);
		c2.resize(nx + 1);
		c3.resize(nx + 1);
		c4.resize(nx + 1);
		c5.resize(nx + 1);

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
			{
				for (int i = 1; i<nx + 1; i++)
				{
					c1[i] = py1[k][j][i];
					c2[i] = py2[k][j][i];
					c3[i] = py3[k][j][i];
					c4[i] = py4[k][j][i];
					c5[i] = py5[k][j][i];
				}
				tdma(ax, bx, ax, c1, c1, nx);
				tdma(ax, bx, ax, c2, c2, nx);
				tdma(ax, bx, ax, c3, c3, nx);
				tdma(ax, bx, ax, c4, c4, nx);
				tdma(ax, bx, ax, c5, c5, nx);
				for (int i = 1; i<nx + 1; i++)
				{
					py1[k][j][i] = c1[i];
					py2[k][j][i] = c2[i];
					py3[k][j][i] = c3[i];
					py4[k][j][i] = c4[i];
					py5[k][j][i] = c5[i];
				}
			}
		c1.resize(ny + 1);
		c2.resize(ny + 1);
		c3.resize(ny + 1);
		c4.resize(ny + 1);
		c5.resize(ny + 1);
		for (int k = 1; k<nz + 1; k++)
			for (int i = 1; i<nx + 1; i++)
			{
				for (int j = 1; j<ny + 1; j++)
				{
					c1[j] = py1[k][j][i];
					c2[j] = py2[k][j][i];
					c3[j] = py3[k][j][i];
					c4[j] = py4[k][j][i];
					c5[j] = py5[k][j][i];
				}
				tdma(ay, by, ay, c1, c1, ny);
				tdma(ay, by, ay, c2, c2, ny);
				tdma(ay, by, ay, c3, c3, ny);
				tdma(ay, by, ay, c4, c4, ny);
				tdma(ay, by, ay, c5, c5, ny);
				for (int j = 1; j<ny + 1; j++)
				{
					py1[k][j][i] = c1[j];
					py2[k][j][i] = c2[j];
					py3[k][j][i] = c3[j];
					py4[k][j][i] = c4[j];
					py5[k][j][i] = c5[j];
				}
			}
		c1.resize(nz + 1);
		c2.resize(nz + 1);
		c3.resize(nz + 1);
		c4.resize(nz + 1);
		c5.resize(nz + 1);
		for (int j = 1; j<ny + 1; j++)
			for (int i = 1; i<nx + 1; i++)
			{
				for (int k = 1; k<nz + 1; k++)
				{
					c1[k] = py1[k][j][i];
					c2[k] = py2[k][j][i];
					c3[k] = py3[k][j][i];
					c4[k] = py4[k][j][i];
					c5[k] = py5[k][j][i];
				}
				tdma(az, bz, az, c1, c1, nz);
				tdma(az, bz, az, c2, c2, nz);
				tdma(az, bz, az, c3, c3, nz);
				tdma(az, bz, az, c4, c4, nz);
				tdma(az, bz, az, c5, c5, nz);
				for (int k = 1; k<nz + 1; k++)
				{
					py1[k][j][i] = c1[k];
					py2[k][j][i] = c2[k];
					py3[k][j][i] = c3[k];
					py4[k][j][i] = c4[k];
					py5[k][j][i] = c5[k];
				}
			}
	}

	void avesa()
	{
		vector<double> ax(nx + 1), bx(nx + 1), ay(ny + 1), by(ny + 1), az(nz + 1), bz(nz + 1), c6;

		for (int i = 0; i<nx + 1; i++)
		{
			ax[i] = -ta;
			bx[i] = 1.0 + 2.0 * ta;
		}

		for (int i = 0; i<ny + 1; i++)
		{
			ay[i] = -ta;
			by[i] = 1.0 + 2.0 * ta;
		}

		for (int i = 0; i<nz + 1; i++)
		{
			az[i] = -ta;
			bz[i] = 1.0 + 2.0 * ta;
		}

		c6.resize(nx + 1);

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
			{
				for (int i = 1; i<nx + 1; i++)
				{
					c6[i] = py6[k][j][i];
				}
				tdma(ax, bx, ax, c6, c6, nx);
				for (int i = 1; i<nx + 1; i++)
				{
					py6[k][j][i] = c6[i];
				}
			}

		c6.resize(ny + 1);
		for (int k = 1; k<nz + 1; k++)
			for (int i = 1; i<nx + 1; i++)
			{
				for (int j = 1; j<ny + 1; j++)
				{
					c6[j] = py6[k][j][i];
				}
				tdma(ay, by, ay, c6, c6, ny);
				for (int j = 1; j<ny + 1; j++)
				{
					py6[k][j][i] = c6[j];
				}
			}

		c6.resize(nz + 1);
		for (int j = 1; j<ny + 1; j++)
			for (int i = 1; i<nx + 1; i++)
			{
				for (int k = 1; k<nz + 1; k++)
				{
					c6[k] = py6[k][j][i];
				}
				tdma(az, bz, az, c6, c6, nz);
				for (int k = 1; k<nz + 1; k++)
				{
					py6[k][j][i] = c6[k];
				}
			}
	}
	//每步R-K推进方法
	void pred(int ims)
	{
		double dtime;

		rrr();

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					py1[k][j][i] = 0;
					py2[k][j][i] = 0;
					py3[k][j][i] = 0;
					py4[k][j][i] = 0;
					py5[k][j][i] = 0;
				}
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					dtime = timl / time[k][j][i];
					py1[k][j][i] = dtime*rr1[k][j][i];
					py2[k][j][i] = dtime*rr2[k][j][i];
					py3[k][j][i] = dtime*rr3[k][j][i];
					py4[k][j][i] = dtime*rr4[k][j][i];
					py5[k][j][i] = dtime*rr5[k][j][i];
					if (isnan(py1[k][j][i])) cout << "py1 error in pred cord" << k << ' ' << j << ' ' << i << ' ' << dtime << ' ' << rr1[k][j][i] << endl;
				}

		if (ims == 1)
			ave();

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					q11[n][k][j][i] = q01[k][j][i] + py1[k][j][i];
					q12[n][k][j][i] = q02[k][j][i] + py2[k][j][i];
					q13[n][k][j][i] = q03[k][j][i] + py3[k][j][i];
					q14[n][k][j][i] = q04[k][j][i] + py4[k][j][i];
					q15[n][k][j][i] = q05[k][j][i] + py5[k][j][i];
					if (isnan(q11[n][k][j][i])) cout << "q11 error in pred cord" << k << ' ' << j << ' ' << i << ' ' << q01[k][j][i] << ' ' << py1[k][j][i] << endl;;
					/*          if(isnan(q12[n][k][j][i])) cout<<"q12 error in pred cord"<<k<<' '<<j<<' '<<i<<' '<<q02[k][j][i]<<' '<<py2[k][j][i]<<endl;
					if(isnan(q13[n][k][j][i])) cout<<"q13 error in pred cord"<<k<<' '<<j<<' '<<i<<' '<<q03[k][j][i]<<' '<<py3[k][j][i]<<endl;
					if(isnan(q14[n][k][j][i])) cout<<"q14 error in pred cord"<<k<<' '<<j<<' '<<i<<' '<<q04[k][j][i]<<' '<<py4[k][j][i]<<endl;
					if(isnan(q15[n][k][j][i])) cout<<"q15 error in pred cord"<<k<<' '<<j<<' '<<i<<' '<<q05[k][j][i]<<' '<<py5[k][j][i]<<endl;*/
				}
	}

	void predsa(int ims)
	{
		double temp = 0;
		double dtime;

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					py6[k][j][i] = 0;
				}

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					dtime = timl / time[k][j][i];
					py6[k][j][i] = dtime*(-qc6[k][j][i] + av6[k][j][i] - ts6[k][j][i] * vvn[layer][k][j][i] + qv6[k][j][i]);
				}

		if (ims == 1)
			avesa();

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					q16[n][k][j][i] = q06[k][j][i] + py6[k][j][i];
				}
	}

	//3步R-K推进
	void march1()
	{
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					q01[k][j][i] = q11[n][k][j][i];
					q02[k][j][i] = q12[n][k][j][i];
					q03[k][j][i] = q13[n][k][j][i];
					q04[k][j][i] = q14[n][k][j][i];
					q05[k][j][i] = q15[n][k][j][i];
				}
		ta = 1.0;
		timl = 0.6*cfl;
		ppp();
		bc();
		step();
		ddd();
		qqq();
		qqqv();
		pred(1);

		timl = 0.6*cfl;
		ppp();
		bc();
		qqq();
		pred(0);

		timl = cfl;
		ppp();
		bc();
		qqq();
		pred(1);
	}
	//粗网格3步R-K推进
	void march2()
	{
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					q01[k][j][i] = q11[n][k][j][i];
					q02[k][j][i] = q12[n][k][j][i];
					q03[k][j][i] = q13[n][k][j][i];
					q04[k][j][i] = q14[n][k][j][i];
					q05[k][j][i] = q15[n][k][j][i];
				}
		ta = 1.0;
		timl = 0.6*cfl;
		ppp();
		bc();
		step();
		dddc();
		qqq();
		qqqv();
		pred(1);

		timl = 0.6*cfl;
		ppp();
		bc();
		qqq();
		pred(0);

		timl = cfl;
		ppp();
		bc();
		qqq();
		pred(1);
	}

	/*
	* ???????????gradscentre(1, pwx, gradc, 1);
	* fotran???????gradscentre(1,pwx(0:nx+1,0:ny+1,0:nz+1),gradc(1,0:nx+1,0:ny+1,0:nz+1))
	*/

	//网格中心梯度
	void gradscentre(int direction, double*** q, double**** dqd, int m)
	{
		double si, flu;

		for (int k = 0; k<nz + 2; k++)
			for (int j = 0; j<ny + 2; j++)
				for (int i = 0; i<nx + 2; i++)
				{
					dqd[m][k][j][i] = 0;
				}

		//********x方向
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 2; i++)
				{
					if (direction == 1)
						si = s2xn[layer][k][j][i];
					else if (direction == 2)
						si = s2yn[layer][k][j][i];
					else if (direction == 3)
						si = s2zn[layer][k][j][i];

					flu = 0.50*(q[k][j][i] + q[k][j][i - 1])*si;
					dqd[m][k][j][i] = dqd[m][k][j][i] + flu;
					dqd[m][k][j][i - 1] = dqd[m][k][j][i - 1] - flu;
				}
		//*********y方向
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 2; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					if (direction == 1)
						si = s3xn[layer][k][j][i];
					else if (direction == 2)
						si = s3yn[layer][k][j][i];
					else if (direction == 3)
						si = s3zn[layer][k][j][i];

					flu = 0.50*(q[k][j - 1][i] + q[k][j][i])*si;
					dqd[m][k][j][i] = dqd[m][k][j][i] + flu;
					dqd[m][k][j - 1][i] = dqd[m][k][j - 1][i] - flu;
				}
		//*********z方向
		for (int k = 1; k<nz + 2; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					if (direction == 1)
						si = s1xn[layer][k][j][i];
					else if (direction == 2)
						si = s1yn[layer][k][j][i];
					else if (direction == 3)
						si = s1zn[layer][k][j][i];

					flu = 0.50*(q[k][j][i] + q[k - 1][j][i])*si;
					dqd[m][k][j][i] = dqd[m][k][j][i] + flu;
					dqd[m][k - 1][j][i] = dqd[m][k - 1][j][i] - flu;
				}

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					dqd[m][k][j][i] = dqd[m][k][j][i] / vvn[layer][k][j][i];
				}
	}

	/*
	* ???????????dsdt(gradc, gradcs, 1)
	* fortran??????? dsdt(gradc(1,0:nx+1,0:ny+1,0:nz+1),gradcs(1,0:nx+1,0:ny+1,0:nz+1))
	*/
	// ????
	void dsdt(double**** q, double**** s, int m)
	{
		double vf, rf, qq1, flu;

		for (int k = 0; k<nz + 2; k++)
			for (int j = 0; j<ny + 2; j++)
				for (int i = 0; i<nx + 2; i++)
				{
					s[m][k][j][i] = 0;
				}

		//******************x方向********************
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 2; i++)
				{
					vx = 0.50*(pvx[k][j][i] + pvx[k][j][i - 1]);
					vy = 0.50*(pvy[k][j][i] + pvy[k][j][i - 1]);
					vz = 0.50*(pvz[k][j][i] + pvz[k][j][i - 1]);
					vf = vx*s2xn[layer][k][j][i] + vy*s2yn[layer][k][j][i] + vz*s2zn[layer][k][j][i];
					rf = rpm*zz2[layer][k][j][i] * s2yn[layer][k][j][i] - rpm*yy2[layer][k][j][i] * s2zn[layer][k][j][i];
					vf = vf + rf;
					qq1 = 0.50*(q[m][k][j][i] + q[m][k][j][i - 1]);
					flu = qq1*vf;
					s[m][k][j][i] = s[m][k][j][i] + flu;
					s[m][k][j][i - 1] = s[m][k][j][i - 1] - flu;
				}

		//******************y方向********************
		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 2; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					vx = 0.50*(pvx[k][j][i] + pvx[k][j - 1][i]);
					vy = 0.50*(pvy[k][j][i] + pvy[k][j - 1][i]);
					vz = 0.50*(pvz[k][j][i] + pvz[k][j - 1][i]);
					vf = vx*s3xn[layer][k][j][i] + vy*s3yn[layer][k][j][i] + vz*s3zn[layer][k][j][i];
					rf = rpm*zz3[layer][k][j][i] * s3yn[layer][k][j][i] - rpm*yy3[layer][k][j][i] * s3zn[layer][k][j][i];
					vf = vf + rf;
					qq1 = 0.50*(q[m][k][j][i] + q[m][k][j - 1][i]);
					flu = qq1*vf;
					s[m][k][j][i] = s[m][k][j][i] + flu;
					s[m][k][j - 1][i] = s[m][k][j - 1][i] - flu;
				}

		//*****************z方向********************
		for (int k = 1; k<nz + 2; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					vx = 0.50*(pvx[k][j][i] + pvx[k - 1][j][i]);
					vy = 0.50*(pvy[k][j][i] + pvy[k - 1][j][i]);
					vz = 0.50*(pvz[k][j][i] + pvz[k - 1][j][i]);
					vf = vx*s1xn[layer][k][j][i] + vy*s1yn[layer][k][j][i] + vz*s1zn[layer][k][j][i];
					rf = rpm*zz1[layer][k][j][i] * s1yn[layer][k][j][i] - rpm*yy1[layer][k][j][i] * s1zn[layer][k][j][i];
					vf = vf + rf;
					qq1 = 0.50*(q[m][k][j][i] + q[m][k - 1][j][i]);
					flu = qq1*vf;
					s[m][k][j][i] = s[m][k][j][i] + flu;
					s[m][k - 1][j][i] = s[m][k - 1][j][i] - flu;
				}

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					s[m][k][j][i] = s[m][k][j][i] / vvn[layer][k][j][i];
				}
	}
	void SAsource()
	{
		// ????
		double*** tur = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		double*** pwx = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		double*** pwy = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		double*** pwz = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		double temp = 0;
		double cvl, tem, fv1, fv2, fv3, rp1, rp2, rp3, w12, w13, w23, w21, w31, w32, ww2, ww, svot, vm, gv;
		double s11, s22, s33, s12, s13, s23, ss2, ss, dd, fr1r1, fr1r2, fr1, yv, ra, ga, fw, dv;

		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
			{
				yy0[k][j][0] = yy0[k][j][1];
				zz0[k][j][0] = zz0[k][j][1];
				yy0[k][j][nx + 1] = yy0[k][j][nx];
				zz0[k][j][nx + 1] = zz0[k][j][nx];
			}

		for (int k = 1; k<nz + 1; k++)
			for (int i = 1; i<nx + 1; i++)
			{
				yy0[k][0][i] = yy0[k][1][i];
				zz0[k][0][i] = zz0[k][1][i];
				yy0[k][ny + 1][i] = yy0[k][ny][i];
				zz0[k][ny + 1][i] = zz0[k][ny][i];
			}

		for (int j = 1; j<ny + 1; j++)
			for (int i = 1; i<nx + 1; i++)
			{
				yy0[0][j][i] = yy0[1][j][i];
				zz0[0][j][i] = zz0[1][j][i];
				yy0[nz + 1][j][i] = yy0[nz][j][i];
				zz0[nz + 1][j][i] = zz0[nz][j][i];
			}

		for (int k = 0; k<nz + 2; k++)
			for (int j = 0; j<ny + 2; j++)
				for (int i = 0; i<nx + 2; i++)
				{
					tur[k][j][i] = q16[n][k][j][i] / q11[n][k][j][i];
					pwx[k][j][i] = pvx[k][j][i];
					pwy[k][j][i] = pvy[k][j][i] + rpm*zz0[k][j][i];
					pwz[k][j][i] = pvz[k][j][i] - rpm*yy0[k][j][i];
				}


		//*******网格中心导数***********
		for (int l = 1; l<13; l++)
			for (int k = 0; k<nz + 2; k++)
				for (int j = 0; j<ny + 2; j++)
					for (int i = 0; i<nx + 2; i++)
					{
						gradc[l][k][j][i] = 0;
					}

		gradscentre(1, pwx, gradc, 1);
		gradscentre(2, pwx, gradc, 2);
		gradscentre(3, pwx, gradc, 3);

		gradscentre(1, pwy, gradc, 4);
		gradscentre(2, pwy, gradc, 5);
		gradscentre(3, pwy, gradc, 6);

		gradscentre(1, pwz, gradc, 7);
		gradscentre(2, pwz, gradc, 8);
		gradscentre(3, pwz, gradc, 9);

		gradscentre(1, tur, gradc, 10);
		gradscentre(2, tur, gradc, 11);
		gradscentre(3, tur, gradc, 12);

		for (int l = 1; l<10; l++)
			for (int k = 0; k<nz + 2; k++)
				for (int j = 0; j<ny + 2; j++)
					for (int i = 0; i<nx + 2; i++)
					{
						gradcs[l][k][j][i] = 0;
					}

		dsdt(gradc, gradcs, 1);
		dsdt(gradc, gradcs, 2);
		dsdt(gradc, gradcs, 3);

		dsdt(gradc, gradcs, 4);
		dsdt(gradc, gradcs, 5);
		dsdt(gradc, gradcs, 6);

		dsdt(gradc, gradcs, 7);
		dsdt(gradc, gradcs, 8);
		dsdt(gradc, gradcs, 9);


		for (int k = 1; k<nz + 1; k++)
			for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					dv = cb2 / sigmav*(pow(gradc[10][k][j][i], 2) + pow(gradc[11][k][j][i], 2) + pow(gradc[12][k][j][i], 2))*q11[n][k][j][i];
					cvl = cvl0*pow(t[k][j][i] / t0, 1.5)*(t0 + ts) / (t[k][j][i] + ts);
					tem = q16[n][k][j][i] / cvl;
					tem = max(tem, pow(10.0, -4));

					fv1 = 1.0 / (1.0 + pow(cv1 / tem, 3));
					fv2 = pow(1.0 + tem / cv2, -3);
					fv3 = (1.0 + tem*fv1)*(1.0 - fv2) / tem;

					w12 = 0.50*(gradc[2][k][j][i] - gradc[4][k][j][i]);
					w13 = 0.50*(gradc[3][k][j][i] - gradc[7][k][j][i]);
					w23 = 0.50*(gradc[6][k][j][i] - gradc[8][k][j][i]);
					ww2 = 4.0*(w12*w12 + w13*w13 + w23*w23);
					ww = sqrt(ww2);

					vm = tur[k][j][i] / (kap*kap*dmini[k][j][i] * dmini[k][j][i]);
					svot = ww*fv3 + vm*fv2;
					gv = cb1*q16[n][k][j][i] * svot;

					s11 = gradc[1][k][j][i];
					s22 = gradc[5][k][j][i];
					s33 = gradc[9][k][j][i];

					s12 = 0.50*(gradc[2][k][j][i] + gradc[4][k][j][i]);
					s13 = 0.50*(gradc[3][k][j][i] + gradc[7][k][j][i]);
					s23 = 0.50*(gradc[6][k][j][i] + gradc[8][k][j][i]);
					ss2 = 4.0*(s12*s12 + s13*s13 + s23*s23) + 2.0*(s11*s11 + s22*s22 + s33*s33);
					ss = sqrt(ss2);

					rp1 = rpm;
					rp2 = 0.0;
					rp3 = 0.0;
					w12 = 0.50*(gradc[2][k][j][i] - gradc[4][k][j][i]) - rp3;
					w13 = 0.50*(gradc[3][k][j][i] - gradc[7][k][j][i]) + rp2;
					w23 = 0.50*(gradc[6][k][j][i] - gradc[8][k][j][i]) - rp1;
					w21 = -w12;
					w31 = -w13;
					w32 = -w23;
					ww2 = 4.0*(w12*w12 + w13*w13 + w23*w23);
					ww = sqrt(ww2);

					fr1r1 = ss / ww;
					dd = max(ss2, 0.09*tur[k][j][i] * tur[k][j][i]);
					fr1r2 = 2.0 / dd / sqrt(dd) / ww*(gradcs[1][k][j][i] * (w12*s12 + w13*s13)
						+ (0.50*(gradcs[2][k][j][i] + gradcs[4][k][j][i]) - s13*rp1)*(w12*s22 + w13*s23 + w21*s11 + w23*s13)
						+ (0.50*(gradcs[3][k][j][i] + gradcs[7][k][j][i]) + s12*rp1)*(w12*s23 + w13*s33 + w31*s11 + w32*s12) + (gradcs[5][k][j][i] - 2.0*s23*rp1)*(w21*s12 + w23*s23)
						+ (0.50*(gradcs[6][k][j][i] + gradcs[8][k][j][i]) + (s22 - s33)*rp1)*(w21*s13 + w23*s33 + w31*s12 + w32*s22) + (gradcs[9][k][j][i] + 2.0*s23*rp1)*(w31*s13 + w32*s23));
					fr1 = (1.0 + cr1)*2.0*fr1r1 / (1.0 + fr1r1)*(1.0 - cr3*atan(cr2*fr1r2)) - cr1;

					gv = gv*fr1;
					ra = vm / svot;
					ga = ra + cw2*(pow(ra, 6) - ra);
					fw = pow((pow(ga, -6) + pow(cw3, -6)) / (1.0 + pow(cw3, -6)), (-1.0 / 6.0));
					yv = cw1*fw*q16[n][k][j][i] * vm*kap*kap;
					qv6[k][j][i] = qv6[k][j][i] + (gv - yv + dv)*vvn[layer][k][j][i];
				}
		free(tur[0][0]);
		free(tur[0]);
		free(tur);
		free(pwx[0][0]);
		free(pwx[0]);
		free(pwx);
		free(pwy[0][0]);
		free(pwy[0]);
		free(pwy);
		free(pwz[0][0]);
		free(pwz[0]);
		free(pwz);
	}

	void qqqvsa()
	{
		double flu6, q16av, tav, cvl;
		double*** tur = malloc3D<double>(nx + 2, ny + 2, nz + 2);


		for (int i = 0; i<nz + 2; i++)
			for (int j = 0; j<ny + 2; j++)
				for (int k = 0; k<nx + 2; k++)
				{
					tur[i][j][k] = q16[n][i][j][k] / q11[n][i][j][k];
				}

		for (int i = 1; i<16; i++)
			for (int j = 0; j<nz + 2; j++)
				for (int k = 0; k<ny + 2; k++)
					for (int l = 1; l<nx + 2; l++)
					{
						gradfi[i][j][k][l] = 0;
					}

		for (int i = 1; i<16; i++)
			for (int j = 0; j<nz + 2; j++)
				for (int k = 1; k<ny + 2; k++)
					for (int l = 0; l<nx + 2; l++)
					{
						gradfj[i][j][k][l] = 0;
					}

		for (int i = 1; i<16; i++)
			for (int j = 1; j<nz + 2; j++)
				for (int k = 0; k<ny + 2; k++)
					for (int l = 0; l<nx + 2; l++)
					{
						gradfk[i][j][k][l] = 0;
					}

		gradsfaceI(s2xn[layer], s3xn[layer], s1xn[layer], tur, gradfi, 13);
		gradsfaceI(s2yn[layer], s3yn[layer], s1yn[layer], tur, gradfi, 14);
		gradsfaceI(s2zn[layer], s3zn[layer], s1zn[layer], tur, gradfi, 15);


		gradsfaceJ(s2xn[layer], s3xn[layer], s1xn[layer], tur, gradfj, 13);
		gradsfaceJ(s2yn[layer], s3yn[layer], s1yn[layer], tur, gradfj, 14);
		gradsfaceJ(s2zn[layer], s3zn[layer], s1zn[layer], tur, gradfj, 15);


		gradsfaceK(s2xn[layer], s3xn[layer], s1xn[layer], tur, gradfk, 13);
		gradsfaceK(s2yn[layer], s3yn[layer], s1yn[layer], tur, gradfk, 14);
		gradsfaceK(s2zn[layer], s3zn[layer], s1zn[layer], tur, gradfk, 15);


		for (int i = 0; i<nz + 2; i++)
			for (int j = 0; j<ny + 2; j++)
				for (int k = 0; k<nx + 2; k++)
				{
					qv6[i][j][k] = 0;
				}


		for (int i = 1; i<nz + 1; i++)
			for (int j = 1; j<ny + 1; j++)
				for (int k = 1; k<nx + 2; k++)
				{
					tav = 0.50*(t[i][j][k - 1] + t[i][j][k]);
					cvl = cvl0*pow((tav / t0), 1.5)*(t0 + ts) / (tav + ts);
					q16av = 0.50*(q16[n][i][j][k - 1] + q16[n][i][j][k]);
					flu6 = (cvl + q16av)*(gradfi[13][i][j][k] * s2xn[layer][i][j][k] + gradfi[14][i][j][k] * s2yn[layer][i][j][k] + gradfi[15][i][j][k] * s2zn[layer][i][j][k]) / sigmav;
					qv6[i][j][k] = qv6[i][j][k] + flu6;
					qv6[i][j][k - 1] = qv6[i][j][k - 1] - flu6;
				}



		for (int i = 1; i<nz + 1; i++)
			for (int j = 1; j<ny + 2; j++)
				for (int k = 1; k<nx + 1; k++)
				{
					tav = 0.50*(t[i][j - 1][k] + t[i][j][k]);
					cvl = cvl0*pow((tav / t0), 1.5)*(t0 + ts) / (tav + ts);
					q16av = 0.50*(q16[n][i][j - 1][k] + q16[n][i][j][k]);
					flu6 = (cvl + q16av)*(gradfj[13][i][j][k] * s3xn[layer][i][j][k] + gradfj[14][i][j][k] * s3yn[layer][i][j][k] + gradfj[15][i][j][k] * s3zn[layer][i][j][k]) / sigmav;
					qv6[i][j][k] = qv6[i][j][k] + flu6;
					qv6[i][j - 1][k] = qv6[i][j - 1][k] - flu6;
				}

		for (int i = 1; i<nz + 2; i++)
			for (int j = 1; j<ny + 1; j++)
				for (int k = 1; k<nx + 1; k++)
				{
					tav = 0.50*(t[i - 1][j][k] + t[i][j][k]);
					cvl = cvl0*pow((tav / t0), 1.5)*(t0 + ts) / (tav + ts);
					q16av = 0.50*(q16[n][i - 1][j][k] + q16[n][i][j][k]);
					flu6 = (cvl + q16av)*(gradfk[13][i][j][k] * s1xn[layer][i][j][k] + gradfk[14][i][j][k] * s1yn[layer][i][j][k] + gradfk[15][i][j][k] * s1zn[layer][i][j][k]) / sigmav;
					qv6[i][j][k] = qv6[i][j][k] + flu6;
					qv6[i - 1][j][k] = qv6[i - 1][j][k] - flu6;
				}
		SAsource();
		free(tur[0][0]);
		free(tur[0]);
		free(tur);
	}

	//SA 3步R-K推进

	void marchsa()
	{
		tsdsa();

		for (int i = 1; i<nz + 1; i++)
			for (int j = 1; j<ny + 1; j++)
				for (int k = 1; k<nx + 1; k++)
				{
					q06[i][j][k] = q16[n][i][j][k];
				}

		ta = 1.5;
		timl = 0.6*cfl;

		ppp();
		// if(myid==0) cout<<"ppp 1 in marchsa"<<endl;
		bc();
		// if(myid==0) cout<<"bc 1 in marchsa"<<endl;
		step();
		// if(myid==0) cout<<"step 1 in marchsa"<<endl;
		dddsa();
		// if(myid==0) cout<<"dddsa 1 in marchsa"<<endl;
		qqqsa();
		// if(myid==0) cout<<"qqqsa 1 in marchsa"<<endl;
		qqqvsa();
		// if(myid==0) cout<<"qqqvsa 1 in marchsa"<<endl;
		predsa(1);
		// if(myid==0) cout<<"predsa 1 in marchsa"<<endl;

		timl = 0.6*cfl;

		ppp();
		//if(myid==0) cout<<"ppp 2 in marchsa"<<endl;
		bc();
		//	if(myid==0) cout<<"bc 2 in marchsa"<<endl;
		qqqsa();
		//if(myid==0) cout<<"qqqsa 2 in marchsa"<<endl;
		predsa(0);
		//if(myid==0) cout<<"ppp 2 in marchsa"<<endl;


		timl = cfl;

		ppp();
		//if(myid==0) cout<<"ppp 3 in marchsa"<<endl;
		bc();
		//if(myid==0) cout<<"bc 3 in marchsa"<<endl;
		qqqsa();
		//if(myid==0) cout<<"qqqsa 1 in marchsa"<<endl;
		predsa(1);
		//if(myid==0) cout<<"predsa 1 in marchsa"<<endl;

		for (int i = 1; i<nz + 1; i++)
			for (int j = 1; j<ny + 1; j++)
				for (int k = 1; k<nx + 1; k++)
				{
					q06[i][j][k] = q16[n][i][j][k];
				}

		timl = 0.125*cfl;

		ppp();
		//if(myid==0) cout<<"ppp 4 in marchsa"<<endl;
		bc();
		//if(myid==0) cout<<"bc 4 in marchsa"<<endl;
		step();
		//if(myid==0) cout<<"step 4 in marchsa"<<endl;
		dddsa();
		//if(myid==0) cout<<"dddsa 4 in marchsa"<<endl;
		qqqsa();
		//if(myid==0) cout<<"qqqsa 4 in marchsa"<<endl;
		qqqvsa();
		//if(myid==0) cout<<"qqqvsa 4 in marchsa"<<endl;
		predsa(0);
		//if(myid==0) cout<<"pred 4 in marchsa"<<endl;
	}

	void residual()
	{
		rms[n] = 0;

		for (int i = 1; i<nz + 1; i++)
			for (int j = 1; j<ny + 1; j++)
				for (int k = 1; k<nx + 1; k++)

				{

					rms[n] = rms[n] + pow((q11[n][i][j][k] - q01[i][j][k])*time[i][j][k] / vvn[layer][i][j][k] / cfl, 2);
					if (isnan(rms[n])) cout << "error in " << i << ' ' << j << ' ' << k << ' ' << q11[n][i][j][k] << ' ' << q01[i][j][k] << ' ' << time[i][j][k] << ' ' << vvn[layer][i][j][k] << ' ' << cfl << endl;

				}

		rms[n] = 0.5*log10(rms[n] / double(nx*ny*nz));
		if (myid == 0) cout << "nitt" << ' ' << nitt << ' ' << "rms" << rms[n] << endl;
		if (rms[n]>rmsm)
			rmsm = rms[n];
	}

	void init()
	{

		int i1, j1, k1;
		double val1, val2, val3, val4, val5, val6, val7, val8, vall;

		for (int i = 1; i<nz + 1; i++)
		{
			i1 = 2 * i;

			for (int j = 1; j<ny + 1; j++)
			{
				j1 = 2 * j;

				for (int k = 1; k<nx + 1; k++)
				{
					k1 = 2 * k;

					val1 = vvn[ign + 1][i1][j1][k1];
					val2 = vvn[ign + 1][i1][j1][k1 - 1];
					val3 = vvn[ign + 1][i1 - 1][j1][k1];
					val4 = vvn[ign + 1][i1 - 1][j1][k1 - 1];
					val5 = vvn[ign + 1][i1][j1 - 1][k1];
					val6 = vvn[ign + 1][i1][j1 - 1][k1 - 1];
					val7 = vvn[ign + 1][i1 - 1][j1 - 1][k1];
					val8 = vvn[ign + 1][i1 - 1][j1 - 1][k1 - 1];
					vall = val1 + val2 + val3 + val4 + val5 + val6 + val7 + val8;

					for (int m = 1; m<nt + 1; m++)
					{
						q11[m][i][j][k] = (q11[m][i1][j1][k1] * val1 + q11[m][i1][j1][k1 - 1] * val2 + q11[m][i1 - 1][j1][k1] * val3 + q11[m][i1 - 1][j1][k1 - 1] * val4 +
							q11[m][i1][j1 - 1][k1] * val5 + q11[m][i1][j1 - 1][k1 - 1] * val6 + q11[m][i1 - 1][j1 - 1][k1] * val7 + q11[m][i1 - 1][j1 - 1][k1 - 1] * val8) / vall;
						q12[m][i][j][k] = (q12[m][i1][j1][k1] * val1 + q12[m][i1][j1][k1 - 1] * val2 + q12[m][i1 - 1][j1][k1] * val3 + q12[m][i1 - 1][j1][k1 - 1] * val4 +
							q12[m][i1][j1 - 1][k1] * val5 + q12[m][i1][j1 - 1][k1 - 1] * val6 + q12[m][i1 - 1][j1 - 1][k1] * val7 + q12[m][i1 - 1][j1 - 1][k1 - 1] * val8) / vall;
						q13[m][i][j][k] = (q13[m][i1][j1][k1] * val1 + q13[m][i1][j1][k1 - 1] * val2 + q13[m][i1 - 1][j1][k1] * val3 + q13[m][i1 - 1][j1][k1 - 1] * val4 +
							q13[m][i1][j1 - 1][k1] * val5 + q13[m][i1][j1 - 1][k1 - 1] * val6 + q13[m][i1 - 1][j1 - 1][k1] * val7 + q13[m][i1 - 1][j1 - 1][k1 - 1] * val8) / vall;
						q14[m][i][j][k] = (q14[m][i1][j1][k1] * val1 + q14[m][i1][j1][k1 - 1] * val2 + q14[m][i1 - 1][j1][k1] * val3 + q14[m][i1 - 1][j1][k1 - 1] * val4 +
							q14[m][i1][j1 - 1][k1] * val5 + q14[m][i1][j1 - 1][k1 - 1] * val6 + q14[m][i1 - 1][j1 - 1][k1] * val7 + q14[m][i1 - 1][j1 - 1][k1 - 1] * val8) / vall;
						q15[m][i][j][k] = (q15[m][i1][j1][k1] * val1 + q15[m][i1][j1][k1 - 1] * val2 + q15[m][i1 - 1][j1][k1] * val3 + q15[m][i1 - 1][j1][k1 - 1] * val4 +
							q15[m][i1][j1 - 1][k1] * val5 + q15[m][i1][j1 - 1][k1 - 1] * val6 + q15[m][i1 - 1][j1 - 1][k1] * val7 + q15[m][i1 - 1][j1 - 1][k1 - 1] * val8) / vall;
						q16[m][i][j][k] = (q16[m][i1][j1][k1] * val1 + q16[m][i1][j1][k1 - 1] * val2 + q16[m][i1 - 1][j1][k1] * val3 + q16[m][i1 - 1][j1][k1 - 1] * val4 +
							q16[m][i1][j1 - 1][k1] * val5 + q16[m][i1][j1 - 1][k1 - 1] * val6 + q16[m][i1 - 1][j1 - 1][k1] * val7 + q16[m][i1 - 1][j1 - 1][k1 - 1] * val8) / vall;
					}



					rr1[i][j][k] = rr1[i1][j1][k1] + rr1[i1][j1][k1 - 1] + rr1[i1 - 1][j1][k1] + rr1[i1 - 1][j1][k1 - 1] +
						rr1[i1][j1 - 1][k1] + rr1[i1 - 1][j1 - 1][k1] + rr1[i1][j1 - 1][k1 - 1] + rr1[i1 - 1][j1 - 1][k1 - 1];
					rr2[i][j][k] = rr2[i1][j1][k1] + rr2[i1 - 1][j1][k1] + rr2[i1][j1][k1 - 1] + rr2[i1 - 1][j1][k1 - 1] +
						rr2[i1][j1 - 1][k1] + rr2[i1 - 1][j1 - 1][k1] + rr2[i1][j1 - 1][k1 - 1] + rr2[i1 - 1][j1 - 1][k1 - 1];
					rr3[i][j][k] = rr3[i1][j1][k1] + rr3[i1 - 1][j1][k1] + rr3[i1][j1][k1 - 1] + rr3[i1 - 1][j1][k1 - 1] +
						rr3[i1][j1 - 1][k1] + rr3[i1 - 1][j1 - 1][k1] + rr3[i1][j1 - 1][k1 - 1] + rr3[i1 - 1][j1 - 1][k1 - 1];
					rr4[i][j][k] = rr4[i1][j1][k1] + rr4[i1 - 1][j1][k1] + rr4[i1][j1][k1 - 1] + rr4[i1 - 1][j1][k1 - 1] +
						rr4[i1][j1 - 1][k1] + rr4[i1 - 1][j1 - 1][k1] + rr4[i1][j1 - 1][k1 - 1] + rr4[i1 - 1][j1 - 1][k1 - 1];
					rr5[i][j][k] = rr5[i1][j1][k1] + rr5[i1 - 1][j1][k1] + rr5[i1][j1][k1 - 1] + rr5[i1 - 1][j1][k1 - 1] +
						rr5[i1][j1 - 1][k1] + rr5[i1 - 1][j1 - 1][k1] + rr5[i1][j1 - 1][k1 - 1] + rr5[i1 - 1][j1 - 1][k1 - 1];
				}
			}
		}
	}

	//粗网格驱动项

	void copr()
	{
		tsd();
		ppp();
		bc();
		step();
		dddc();
		qqq();
		qqqv();

		for (int i = 1; i<nz + 1; i++)
			for (int j = 1; j<ny + 1; j++)
				for (int k = 1; k<nx + 1; k++)
				{
					qp1[i][j][k] = rr1[i][j][k] + qc1[i][j][k] - av1[i][j][k] + ts1[i][j][k] * vvn[layer][i][j][k];
					qp2[i][j][k] = rr2[i][j][k] + qc2[i][j][k] - av2[i][j][k] + ts2[i][j][k] * vvn[layer][i][j][k] - qv2[i][j][k];
					qp3[i][j][k] = rr3[i][j][k] + qc3[i][j][k] - av3[i][j][k] + ts3[i][j][k] * vvn[layer][i][j][k] - qv3[i][j][k];
					qp4[i][j][k] = rr4[i][j][k] + qc4[i][j][k] - av4[i][j][k] + ts4[i][j][k] * vvn[layer][i][j][k] - qv4[i][j][k];
					qp5[i][j][k] = rr5[i][j][k] + qc5[i][j][k] - av5[i][j][k] + ts5[i][j][k] * vvn[layer][i][j][k] - qv5[i][j][k];
				}
	}

	double dp(double a, double b, double c)
	{
		return (3.0*a - 2.0*b - c) / 64.0;
	}

	//插值

	void update()
	{
		int ig, i1, j1, k1, iw, ie, jn, js, kf, kb;
		double dpy, dpiw, dpie, dpjn, dpjs, dpkf, dpkb;

		for (ig = ign; ig<nng; ig++)
		{
			nx = nnx[ig];
			ny = nny[ig];
			nz = nnz[ig];

			for (int k = 1; k<nz + 1; k++)
			{
				/*		i1 = 2 * i - 1;
				iw = max(1, i - 1);
				ie = min(nx, i + 1);*/

				k1 = 2 * k - 1;
				kf = max(1, k - 1);
				kb = min(nz, k + 1);

				for (int j = 1; j<ny + 1; j++)
				{
					j1 = 2 * j - 1;
					js = max(1, j - 1);
					jn = min(ny, j + 1);

					for (int i = 1; i<nx + 1; i++)
					{
						/*		k1 = 2 * k - 1;
						kf = max(1, k - 1);
						kb = min(nz, k + 1);*/

						i1 = 2 * i - 1;
						iw = max(1, i - 1);
						ie = min(nx, i + 1);



						dpy = py1[k][j][i];
						dpiw = dp(py1[k][j][iw], dpy, py1[k][j][ie]);
						dpie = dp(py1[k][j][ie], dpy, py1[k][j][iw]);
						dpjs = dp(py1[k][js][i], dpy, py1[k][jn][i]);
						dpjn = dp(py1[k][jn][i], dpy, py1[k][js][i]);

						dpkf = dp(py1[kf][j][i], dpy, py1[kb][j][i]);
						dpkb = dp(py1[kb][j][i], dpy, py1[kf][j][i]);

						q01[k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
						q01[k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
						q01[k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
						q01[k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
						q01[k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
						q01[k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
						q01[k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
						q01[k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;



						dpy = py2[k][j][i];

						dpiw = dp(py2[k][j][iw], dpy, py2[k][j][ie]);
						dpie = dp(py2[k][j][ie], dpy, py2[k][j][iw]);
						dpjs = dp(py2[k][js][i], dpy, py2[k][jn][i]);
						dpjn = dp(py2[k][jn][i], dpy, py2[k][js][i]);

						dpkf = dp(py2[kf][j][i], dpy, py2[kb][j][i]);
						dpkb = dp(py2[kb][j][i], dpy, py2[kf][j][i]);

						q02[k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
						q02[k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
						q02[k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
						q02[k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
						q02[k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
						q02[k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
						q02[k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
						q02[k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;



						dpy = py3[k][j][i];

						dpiw = dp(py3[k][j][iw], dpy, py3[k][j][ie]);
						dpie = dp(py3[k][j][ie], dpy, py3[k][j][iw]);
						dpjs = dp(py3[k][js][i], dpy, py3[k][jn][i]);
						dpjn = dp(py3[k][jn][i], dpy, py3[k][js][i]);

						dpkf = dp(py3[kf][j][i], dpy, py3[kb][j][i]);
						dpkb = dp(py3[kb][j][i], dpy, py3[kf][j][i]);

						q03[k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
						q03[k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
						q03[k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
						q03[k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
						q03[k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
						q03[k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
						q03[k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
						q03[k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;



						dpy = py4[k][j][i];

						dpiw = dp(py4[k][j][iw], dpy, py4[k][j][ie]);
						dpie = dp(py4[k][j][ie], dpy, py4[k][j][iw]);
						dpjs = dp(py4[k][js][i], dpy, py4[k][jn][i]);
						dpjn = dp(py4[k][jn][i], dpy, py4[k][js][i]);

						dpkf = dp(py4[kf][j][i], dpy, py4[kb][j][i]);
						dpkb = dp(py4[kb][j][i], dpy, py4[kf][j][i]);

						q04[k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
						q04[k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
						q04[k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
						q04[k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
						q04[k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
						q04[k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
						q04[k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
						q04[k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;



						dpy = py5[k][j][i];

						dpiw = dp(py5[k][j][iw], dpy, py5[k][j][ie]);
						dpie = dp(py5[k][j][ie], dpy, py5[k][j][iw]);
						dpjs = dp(py5[k][js][i], dpy, py5[k][jn][i]);
						dpjn = dp(py5[k][jn][i], dpy, py5[k][js][i]);

						dpkf = dp(py5[kf][j][i], dpy, py5[kb][j][i]);
						dpkb = dp(py5[kb][j][i], dpy, py5[kf][j][i]);

						q05[k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
						q05[k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
						q05[k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
						q05[k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
						q05[k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
						q05[k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
						q05[k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
						q05[k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;
					}
				}
			}

			nx = nnx[ig + 1];
			ny = nny[ig + 1];
			nz = nnz[ig + 1];

			/*			py1=q01;
			py2=q02;
			py3=q03;
			py4=q04;
			py5=q05;
			*/
			for (int i = 1; i<nz + 1; i++)
				for (int j = 1; j<ny + 1; j++)
					for (int k = 1; k<nx + 1; k++)
					{
						py1[i][j][k] = q01[i][j][k];
						py2[i][j][k] = q02[i][j][k];
						py3[i][j][k] = q03[i][j][k];
						py4[i][j][k] = q04[i][j][k];
						py5[i][j][k] = q05[i][j][k];
					}
		}

		for (int i = 1; i<nz + 1; i++)
			for (int j = 1; j<ny + 1; j++)
				for (int k = 1; k<nx + 1; k++)
				{
					q31[n][i][j][k] = q31[n][i][j][k] + py1[i][j][k];
					q32[n][i][j][k] = q32[n][i][j][k] + py2[i][j][k];
					q33[n][i][j][k] = q33[n][i][j][k] + py3[i][j][k];
					q34[n][i][j][k] = q34[n][i][j][k] + py4[i][j][k];
					q35[n][i][j][k] = q35[n][i][j][k] + py5[i][j][k];
				}

		nx = nnx[ign];
		ny = nny[ign];
		nz = nnz[ign];
	}

	//插值跳阶段

	void update2()
	{
		int i1, j1, k1, iw, ie, jn, js, kf, kb;
		double dpy, dpiw, dpie, dpjn, dpjs, dpkf, dpkb;

		nx = nnx[nng];
		ny = nny[nng];
		nz = nnz[nng];

		for (n = 1; n<nt + 1; n++)
		{
			for (int i = 1; i<nz + 1; i++)
				for (int j = 1; j<ny + 1; j++)
					for (int k = 1; k<nx + 1; k++)
					{
						py1[i][j][k] = q11[n][i][j][k];
						py2[i][j][k] = q12[n][i][j][k];
						py3[i][j][k] = q13[n][i][j][k];
						py4[i][j][k] = q14[n][i][j][k];
						py5[i][j][k] = q15[n][i][j][k];
						py6[i][j][k] = q16[n][i][j][k];
					}

			for (int k = 1; k<nz + 1; k++)
			{
				/*				i1 = 2 * i - 1;
				iw = max(1, i - 1);
				ie = min(nx, i + 1);*/

				k1 = 2 * k - 1;
				kf = max(1, k - 1);
				kb = min(nz, k + 1);

				for (int j = 1; j<ny + 1; j++)
				{
					j1 = 2 * j - 1;
					js = max(1, j - 1);
					jn = min(ny, j + 1);

					for (int i = 1; i<nx + 1; i++)
					{
						/*k1 = 2 * k - 1;
						kf = max(1, k - 1);
						kb = min(nz, k + 1);
						*/
						i1 = 2 * i - 1;
						iw = max(1, i - 1);
						ie = min(nx, i + 1);


						dpy = py1[k][j][i];

						dpiw = dp(py1[k][j][iw], dpy, py1[k][j][ie]);
						dpie = dp(py1[k][j][ie], dpy, py1[k][j][iw]);
						dpjs = dp(py1[k][js][i], dpy, py1[k][jn][i]);
						dpjn = dp(py1[k][jn][i], dpy, py1[k][js][i]);
						dpkf = dp(py1[kf][j][i], dpy, py1[kb][j][i]);
						dpkb = dp(py1[kb][j][i], dpy, py1[kf][j][i]);

						q11[n][k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
						q11[n][k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
						q11[n][k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
						q11[n][k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
						q11[n][k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
						q11[n][k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
						q11[n][k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
						q11[n][k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;



						dpy = py2[k][j][i];
						dpiw = dp(py2[k][j][iw], dpy, py2[k][j][ie]);
						dpie = dp(py2[k][j][ie], dpy, py2[k][j][iw]);
						dpjs = dp(py2[k][js][i], dpy, py2[k][jn][i]);
						dpjn = dp(py2[k][jn][i], dpy, py2[k][js][i]);
						dpkf = dp(py2[kf][j][i], dpy, py2[kb][j][i]);
						dpkb = dp(py2[kb][j][i], dpy, py2[kf][j][i]);

						q12[n][k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
						q12[n][k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
						q12[n][k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
						q12[n][k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
						q12[n][k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
						q12[n][k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
						q12[n][k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
						q12[n][k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;



						dpy = py3[k][j][i];
						dpiw = dp(py3[k][j][iw], dpy, py3[k][j][ie]);
						dpie = dp(py3[k][j][ie], dpy, py3[k][j][iw]);
						dpjs = dp(py3[k][js][i], dpy, py3[k][jn][i]);
						dpjn = dp(py3[k][jn][i], dpy, py3[k][js][i]);
						dpkf = dp(py3[kf][j][i], dpy, py3[kb][j][i]);
						dpkb = dp(py3[kb][j][i], dpy, py3[kf][j][i]);

						q13[n][k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
						q13[n][k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
						q13[n][k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
						q13[n][k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
						q13[n][k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
						q13[n][k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
						q13[n][k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
						q13[n][k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;



						dpy = py4[k][j][i];
						dpiw = dp(py4[k][j][iw], dpy, py4[k][j][ie]);
						dpie = dp(py4[k][j][ie], dpy, py4[k][j][iw]);
						dpjs = dp(py4[k][js][i], dpy, py4[k][jn][i]);
						dpjn = dp(py4[k][jn][i], dpy, py4[k][js][i]);
						dpkf = dp(py4[kf][j][i], dpy, py4[kb][j][i]);
						dpkb = dp(py4[kb][j][i], dpy, py4[kf][j][i]);

						q14[n][k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
						q14[n][k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
						q14[n][k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
						q14[n][k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
						q14[n][k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
						q14[n][k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
						q14[n][k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
						q14[n][k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;



						dpy = py5[k][j][i];
						dpiw = dp(py5[k][j][iw], dpy, py5[k][j][ie]);
						dpie = dp(py5[k][j][ie], dpy, py5[k][j][iw]);
						dpjs = dp(py5[k][js][i], dpy, py5[k][jn][i]);
						dpjn = dp(py5[k][jn][i], dpy, py5[k][js][i]);
						dpkf = dp(py5[kf][j][i], dpy, py5[kb][j][i]);
						dpkb = dp(py5[kb][j][i], dpy, py5[kf][j][i]);

						q15[n][k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
						q15[n][k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
						q15[n][k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
						q15[n][k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
						q15[n][k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
						q15[n][k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
						q15[n][k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
						q15[n][k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;



						dpy = py6[k][j][i];
						dpiw = dp(py6[k][j][iw], dpy, py6[k][j][ie]);
						dpie = dp(py6[k][j][ie], dpy, py6[k][j][iw]);
						dpjs = dp(py6[k][js][i], dpy, py6[k][jn][i]);
						dpjn = dp(py6[k][jn][i], dpy, py6[k][js][i]);
						dpkf = dp(py6[kf][j][i], dpy, py6[kb][j][i]);
						dpkb = dp(py6[kb][j][i], dpy, py6[kf][j][i]);

						q16[n][k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
						q16[n][k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
						q16[n][k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
						q16[n][k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
						q16[n][k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
						q16[n][k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
						q16[n][k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
						q16[n][k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;
					}
				}
			}
		}
	}

	//相对马赫数

	void wmma()
	{
		double qqw, a;

		for (int i = 1; i<nz + 1; i++)
			for (int j = 1; j<ny + 1; j++)
				for (int k = 1; k<nx + 1; k++)
				{
					y1 = yy0[i][j][k];
					z1 = zz0[i][j][k];

					vx = pvx[i][j][k];
					vy = pvy[i][j][k];
					vz = pvz[i][j][k];

					wx = vx;
					wy = vy + rpm*z1;
					wz = vz - rpm*y1;

					qqw = wx*wx + wy*wy + wz*wz;
					a = 1.40*p[i][j][k] / q11[n][i][j][k];

					wma[i][j][k] = sqrt(qqw / a);
				}
	}

	//进出口流量

	void flow()
	{
		double  vf1, vf2, qin, qout, qinn, qoutt;
		vector<double> qin0, qout0;

		qin = 0;
		qout = 0;

		for (n = 1; n<nt + 1; n++)
			for (int j = 1; j<nz + 1; j++)
				for (int k = 1; k<ny + 1; k++)
					//	for (n = 1; n<nt + 1; n++)
				{
					vf1 = -(s2xn[nng][j][k][1] * q12[n][j][k][1] + s2yn[nng][j][k][1] * q13[n][j][k][1] + s2zn[nng][j][k][1] * q14[n][j][k][1]);
					qin = qin + vf1;

					vf2 = -(s2xn[nng][j][k][nx + 1] * q12[n][j][k][nx] + s2yn[nng][j][k][nx + 1] * q13[n][j][k][nx] + s2zn[nng][j][k][nx + 1] * q14[n][j][k][nx]);
					qout = qout + vf2;
				}

		qin = qin / double(nt);
		qout = qout / double(nt);

		if (myid == 0)
		{

			qin0.resize(lbb);
			qout0.resize(lbb);

			qin0[0] = qin;
			qinn = 0;

			qout0[0] = qout;
			qoutt = 0;

			for (int j = 1; j<lbb; j++)
			{
				MPI_Recv(&qin0[j], 1, MPI_DOUBLE, j, 3, MPI_COMM_WORLD, &status);
				MPI_Recv(&qout0[j], 1, MPI_DOUBLE, j, 4, MPI_COMM_WORLD, &status);
			}

			for (int j = 0; j<lbb; j++)
			{
				qinn = qinn + qin0[j];
				qoutt = qoutt + qout0[j];
			}

			f1 << setw(5) << nitt;
			f1 << setw(15) << setprecision(6) << qinn;
			f2 << setw(5) << nitt;
			f2 << setw(15) << setprecision(6) << qoutt;
		}
		else
		{
			MPI_Send(&qin, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
			MPI_Send(&qout, 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
		}
	}

	//进出口量

	void inlout()
	{
		int mml;
		double  vf1, vf2, qin, qout1, qq2, h0, t2, p2, pout, tout, d1, d2, eff, temp, ftt, fpp;
		vector<double>rr0, qout, fp, ft, qin0, qout0, d10, d20, eff0;
		double qinn, qoutt, d11, d22, efff;
		fstream f4, f5, f8, f9, f11;
		string zonename;

		if (myid == 0)
		{
			f4.open("inflow.dat", std::fstream::in | std::fstream::out | std::fstream::app);
			f5.open("outflow.dat", std::fstream::in | std::fstream::out | std::fstream::app);
			f8.open("pbi.dat", std::fstream::in | std::fstream::out | std::fstream::app);
			f9.open("tbi.dat", std::fstream::in | std::fstream::out | std::fstream::app);
			f11.open("eff.dat", std::fstream::in | std::fstream::out | std::fstream::app);
		}

		for (n = 1; n<nt + 1; n++)
		{
			qin = 0;

			for (int j = 1; j<nz + 1; j++)
				for (int k = 1; k<ny + 1; k++)
				{
					vf1 = -(s2xn[layer][j][k][1] * q12[n][j][k][1] + s2yn[layer][j][k][1] * q13[n][j][k][1] + s2zn[layer][j][k][1] * q14[n][j][k][1]);
					qin = qin + vf1;
				}

			if (myid == 0)
			{
				qin0.resize(lbb);
				qin0[0] = qin;
				qinn = 0;

				for (int j = 1; j<lbb; j++)
					MPI_Recv(&qin0[j], 1, MPI_DOUBLE, j, 3, MPI_COMM_WORLD, &status);

				for (int j = 0; j<lbb; j++)
					qinn = qinn + qin0[j];

				f4 << setw(6) << setprecision(17) << double(n - 1) / double(nt);
				f4 << setw(15) << setprecision(17) << qinn;
			}
			else
				MPI_Send(&qin, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
		}

		rr0.resize(ny + 1);
		qout.resize(ny + 1);
		fp.resize(ny + 1);
		ft.resize(ny + 1);

		for (int j = 1; j<ny + 1; j++)
		{
			y1 = yy[layer][1][j][nx];
			z1 = zz[layer][1][j][nx];

			rr0[j] = sqrt(y1*y1 + z1*z1);
		}

		for (n = 1; n<nt + 1; n++)
		{
			for (int j = 1; j<ny + 1; j++)
			{
				qout[j] = 0.0;
				fp[j] = 0.0;
				ft[j] = 0.0;
			}

			qout1 = 0.0;
			fpp = 0.0;
			ftt = 0.0;

			ostringstream os;
			os << n;
			zonename = os.str();


			fstream f44, f55, f66;
			string file = "pbispan-timel" + trim(zonename) + ".dat";
			f44.open(file.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

			file = "tbispan-timel" + trim(zonename) + ".dat";
			f55.open(file.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

			file = "effspan-timel" + trim(zonename) + ".dat";
			f66.open(file.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

			for (int j = 1; j<ny + 1; j++)
			{
				temp = (rr0[j] - rr0[1]) / (rr0[ny] - rr0[1])*100.0;

				for (int k = 1; k<nz + 1; k++)
				{
					vf2 = -(s2xn[layer][k][j][nx + 1] * q12[n][k][j][nx] + s2yn[layer][k][j][nx + 1] * q13[n][k][j][nx] + s2zn[layer][k][j][nx + 1] * q14[n][k][j][nx]);
					vx = q12[n][k][j][nx] / q11[n][k][j][nx];
					vy = q13[n][k][j][nx] / q11[n][k][j][nx];
					vz = q14[n][k][j][nx] / q11[n][k][j][nx];

					qq2 = vx*vx + vy*vy + vz*vz;
					pp = 0.40*(q15[n][k][j][nx] - 0.50*q11[n][k][j][nx] * qq2);
					h0 = (q15[n][k][j][nx] + pp) / q11[n][k][j][nx];
					t2 = h0 / cp;
					p2 = pp*pow(1.0 - 0.50*qq2 / h0, -3.5);

					ft[j] = ft[j] + t2*vf2;
					fp[j] = fp[j] + p2*vf2;
					qout[j] = qout[j] + vf2;
				}

				pout = fp[j] / qout[j];
				tout = ft[j] / qout[j];

				d1 = pout / pt;
				d2 = tout / (ht / cp);
				eff = (pow(d1, 2.0 / 7.0) - 1.0) / (d2 - 1.0)*100.0;

				if (myid == 0)
				{
					d10.resize(lbb);
					d20.resize(lbb);
					eff0.resize(lbb);

					d10[0] = d1;
					d20[0] = d2;
					eff0[0] = eff;

					d11 = 0.0;
					d22 = 0.0;
					efff = 0.0;

					for (int i = 1; i<lbb; i++)
					{
						MPI_Recv(&d10[i], 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status);
						MPI_Recv(&d20[i], 1, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, &status);
						MPI_Recv(&eff0[i], 1, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, &status);
					}

					for (int i = 0; i<lbb; i++)
					{
						d11 = d11 + d10[i];
						d22 = d22 + d20[i];
						efff = efff + eff0[i];
					}

					f44 << setw(15) << setprecision(6) << d11 / (mml*lbb) << "         ";
					f44 << setw(15) << setprecision(6) << temp << endl;
					f55 << setw(15) << setprecision(6) << d11 / (mml*lbb) << "        ";
					f55 << setw(15) << setprecision(6) << temp << endl;
					f66 << setw(15) << setprecision(6) << d11 / (mml*lbb) << "        ";
					f66 << setw(15) << setprecision(6) << temp << endl;
				}
				else
				{
					MPI_Send(&d1, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
					MPI_Send(&d2, 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
					MPI_Send(&eff, 1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
				}

				ftt = ftt + ft[j];
				fpp = fpp + fp[j];
				qout1 = qout1 + qout[j];

			}

			f44.close();
			f55.close();
			f66.close();

			pout = fpp / qout1;
			tout = ftt / qout1;

			d1 = pout / pt;
			d2 = tout / (ht / cp);
			eff = (pow(d1, 2.0 / 7.0) - 1.0) / (d2 - 1.0)*100.0;

			if (myid == 0)
			{
				qout0.resize(lbb);
				d10.resize(lbb);
				d20.resize(lbb);
				eff0.resize(lbb);

				qout0[0] = qout1;
				d10[0] = d1;
				d20[0] = d2;
				eff0[0] = eff;

				qoutt = 0.0;
				d11 = 0.0;
				d22 = 0.0;
				efff = 0.0;

				for (int j = 1; j<lbb; j++)
				{

					MPI_Recv(&qout0[j], 1, MPI_DOUBLE, j, 2, MPI_COMM_WORLD, &status);
					MPI_Recv(&d10[j], 1, MPI_DOUBLE, j, 3, MPI_COMM_WORLD, &status);
					MPI_Recv(&d20[j], 1, MPI_DOUBLE, j, 4, MPI_COMM_WORLD, &status);
					MPI_Recv(&eff0[j], 1, MPI_DOUBLE, j, 5, MPI_COMM_WORLD, &status);
				}

				for (int j = 0; j<lbb; j++)
				{
					qoutt = qoutt + qout0[j];
					d11 = d11 + d10[j];
					d22 = d22 + d20[j];
					efff = efff + eff0[j];
				}

				f5 << setw(6) << setprecision(3) << double(n - 1) / double(nt);
				f5 << setw(15) << setprecision(6) << qoutt;
				f8 << setw(6) << setprecision(3) << double(n - 1) / double(nt);
				f8 << setw(15) << setprecision(6) << d11 / lbb;
				f9 << setw(6) << setprecision(3) << double(n - 1) / double(nt);
				f9 << setw(15) << setprecision(6) << d22 / lbb;
				f11 << setw(6) << setprecision(3) << double(n - 1) / double(nt);
				f11 << setw(15) << setprecision(6) << efff / lbb;
			}
			else
			{
				MPI_Send(&qout1, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
				MPI_Send(&d1, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
				MPI_Send(&d2, 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
				MPI_Send(&eff, 1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
			}
		}
	}

	//叶片

	void lbout()
	{
		//open(44,file='out-'//trim(adjustl(id_m))//'myid.dat')  !// ?????????????'abc' // '.txt' ?????'abc.txt',Trim ??????id_m??????adjustl ??????id_m??????
		//write(44,*) "VARIABLES=x,y,z,u,v,w,wma,pressure,density,ut"

		string zonename;
		fstream f44;
		string file = "out-" + trim(id_m) + "myid.dat";

		f44.open(file.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
		f44 << "VARIABLES=x,y,z,u,v,w,wma,pressure,density,ut";

		for (n = 1; n<nt + 1; n++)
		{
			ppp();
			wmma();

			ostringstream os;
			os << n;
			zonename = os.str();
			zonename = "out_" + trim(zonename);

			if (n == 1)
			{
				f44 << "ZONE T=" << zonename << "I=" << nx + 1 << "J=" << ny + 1 << "K=" << nz + 1
					<< "DATAPACKING=BLOCK, VARLOCATION=([4-10]=CELLCENTERED),SOLUTIONTIME="
					<< setw(8) << setprecision(4) << double(n - 1) / double(nt) << endl;

				int count = 0;

				for (int i = 1; i<nz + 2; i++)
					for (int j = 1; j<ny + 2; j++)
						for (int k = 1; k<nx + 2; k++)
						{
							f44 << setprecision(15) << x[i][j][k] << ' ';
							count++;
							if (count == 3)
							{
								f44 << endl;
								count = 0;
							}
						}

				f44 << endl;

				count = 0;
				for (int i = 1; i<nz + 2; i++)
					for (int j = 1; j<ny + 2; j++)
						for (int k = 1; k<nx + 2; k++)
						{
							f44 << setprecision(15) << y[i][j][k] << ' ';
							count++;
							if (count == 3)
							{
								f44 << endl;
								count = 0;
							}
						}

				count = 0;
				f44 << endl;

				for (int i = 1; i<nz + 2; i++)
					for (int j = 1; j<ny + 2; j++)
						for (int k = 1; k<nx + 2; k++)
						{
							f44 << setprecision(15) << z[i][j][k] << ' ';
							count++;
							if (count == 3)
							{
								f44 << endl;
								count = 0;
							}
						}

				count = 0;
				f44 << endl;

				for (int i = 1; i<nz + 1; i++)
					for (int j = 1; j<ny + 1; j++)
						for (int k = 1; k<nx + 1; k++)
						{
							f44 << setprecision(15) << pvx[i][j][k] << ' ';
							count++;
							if (count == 3)
							{
								f44 << endl;
								count = 0;
							}
						}

				count = 0;
				f44 << endl;

				for (int i = 1; i<nz + 1; i++)
					for (int j = 1; j<ny + 1; j++)
						for (int k = 1; k<nx + 1; k++)
						{
							f44 << setprecision(15) << pvy[i][j][k] << ' ';
							count++;
							if (count == 3)
							{
								f44 << endl;
								count = 0;
							}
						}

				count = 0;
				f44 << endl;

				for (int i = 1; i<nz + 1; i++)
					for (int j = 1; j<ny + 1; j++)
						for (int k = 1; k<nx + 1; k++)
						{
							f44 << setprecision(15) << pvz[i][j][k] << ' ';
							count++;
							if (count == 3)
							{
								f44 << endl;
								count = 0;
							}
						}

				count = 0;
				f44 << endl;

				for (int i = 1; i<nz + 1; i++)
					for (int j = 1; j<ny + 1; j++)
						for (int k = 1; k<nx + 1; k++)
						{
							f44 << setprecision(15) << wma[i][j][k] << ' ';
							count++;
							if (count == 3)
							{
								f44 << endl;
								count = 0;
							}
						}

				count = 0;
				f44 << endl;

				for (int i = 1; i<nz + 1; i++)
					for (int j = 1; j<ny + 1; j++)
						for (int k = 1; k<nx + 1; k++)
						{
							f44 << setprecision(15) << p[i][j][k] << ' ';
							count++;
							if (count == 3)
							{
								f44 << endl;
								count = 0;
							}
						}

				f44 << endl;
				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int j = 1; j<ny + 1; j++)
						for (int k = 1; k<nx + 1; k++)
						{
							f44 << setprecision(15) << q11[n][i][j][k] << ' ';
							count++;
							if (count == 3)
							{
								f44 << endl;
								count = 0;
							}
						}

				f44 << endl;
				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int j = 1; j<ny + 1; j++)
						for (int k = 1; k<nx + 1; k++)
						{
							f44 << setprecision(15) << q16[n][i][j][k] << ' ';
							count++;
							if (count == 3)
							{
								f44 << endl;
								count = 0;
							}
						}

				count = 0;
				f44 << endl;
			}
			else
			{
				f44 << "ZONE T=" << zonename << "I=" << nx + 1 << "J=" << ny + 1 << "K=" << nz + 1
					<< "DATAPACKING=BLOCK, VARSHARELIST=([1-3]=1),VARLOCATION=([4-10]=CELLCENTERED),SOLUTIONTIME="
					<< setw(8) << setprecision(4) << double(n - 1) / double(nt);

				int count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int j = 1; j<ny + 1; j++)
						for (int k = 1; k<nx + 1; k++)
						{
							f44 << setprecision(15) << pvx[i][j][k] << ' ';
							count++;
							if (count == 3)
							{
								f44 << endl;
								count = 0;
							}
						}

				count = 0;
				f44 << endl;

				for (int i = 1; i<nz + 1; i++)
					for (int j = 1; j<ny + 1; j++)
						for (int k = 1; k<nx + 1; k++)
						{
							f44 << setprecision(15) << pvy[i][j][k] << ' ';
							count++;
							if (count == 3)
							{
								f44 << endl;
								count = 0;
							}
						}

				count = 0;
				f44 << endl;

				for (int i = 1; i<nz + 1; i++)
					for (int j = 1; j<ny + 1; j++)
						for (int k = 1; k<nx + 1; k++)
						{
							f44 << setprecision(15) << pvz[i][j][k] << ' ';
							count++;
							if (count == 3)
							{
								f44 << endl;
								count = 0;
							}
						}

				f44 << endl;
				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int j = 1; j<ny + 1; j++)
						for (int k = 1; k<nx + 1; k++)
						{
							f44 << setprecision(15) << wma[i][j][k] << ' ';
							count++;
							if (count == 3)
							{
								f44 << endl;
								count = 0;
							}
						}

				f44 << endl;
				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int j = 1; j<ny + 1; j++)
						for (int k = 1; k<nx + 1; k++)
						{
							f44 << setprecision(15) << p[i][j][k] << ' ';
							count++;
							if (count == 3)
							{
								f44 << endl;
								count = 0;
							}
						}

				f44 << endl;
				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int j = 1; j<ny + 1; j++)
						for (int k = 1; k<nx + 1; k++)
						{
							f44 << setprecision(15) << q11[n][i][j][k] << ' ';
							count++;
							if (count == 3)
							{
								f44 << endl;
								count = 0;
							}
						}

				f44 << endl;
				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int j = 1; j<ny + 1; j++)
						for (int k = 1; k<nx + 1; k++)
						{
							f44 << setprecision(15) << q16[n][i][j][k] << ' ';
							count++;
							if (count == 3)
							{
								f44 << endl;
								count = 0;
							}
						}
			}
		}

		f44.close();
	}

	//三次样条插值

	double wl(vector<double>medx, vector<double>medr, double spax, int n1, int n2)
	{
		vector<double> h, g, lamt, miu, a, m;

		h.resize(n1);
		g.resize(n1);
		lamt.resize(n1);
		miu.resize(n1);
		a.resize(n1);
		m.resize(n1 + 1);

		double h1, spar;

		spar = 0;

		m[1] = (medr[2] - medr[1]) / (medx[2] - medx[1]);
		m[n1] = (medr[n1] - medr[n1 - 1]) / (medx[n1] - medx[n1 - 1]);

		for (int i = 1; i<n1; i++)
		{
			h[i] = medx[i + 1] - medx[i];
		}

		for (int i = 2; i<n1; i++)
		{
			a[i] = 2.0;
			lamt[i] = h[i] / (h[i] + h[i - 1]);
			miu[i] = 1 - lamt[i];
			g[i] = 3.0*(lamt[i] * (medr[i] - medr[i - 1]) / h[i - 1] + miu[i] * (medr[i + 1] - medr[i]) / h[i]);
		}

		g[2] = g[2] - lamt[2] * m[1];
		g[n1 - 1] = g[n1 - 1] - miu[n1 - 1] * m[n1];

		vector<double> u(n1 - 1), l(n1 - 1), y(n1 - 1);

		u[2] = a[2];
		y[2] = g[2];

		for (int i = 3; i<n1 - 1; i++)
		{
			l[i] = lamt[i] / u[i - 1];
			u[i] = a[i] - l[i] * miu[i - 1];
			y[i] = g[i] - l[i] * y[i - 1];
		}

		m[n1 - 2] = y[n] / u[n];

		for (int i = n1 - 2; i>0; i--)
		{
			m[i] = (y[i] - miu[i - 1] * m[i + 1]) / u[i];
		}
		//////////////////////////////////////////////////////////////////

		for (int j = 1; j<n2 + 1; j++)
		{
			int i;
			if (spax<medx[n1])
			{
				i = 1;
				while (spax>medx[i + 1])
				{
					i++;
				}
			}
			else
				i = n1 - 1;

			h1 = (spax - medx[i]) / h[i];
			spar = spar + medr[i] * (2.0*h1 + 1.0)*(h1 - 1.0)*(h1 - 1.0) + m[i] * h[i] * h1*(h1 - 1.0)*(h1 - 1.0);
			h1 = (medx[i + 1] - spax) / h[i];
			spar = spar + medr[i + 1] * (2.0*h1 + 1.0)*(h1 - 1.0)*(h1 - 1.0) - m[i + 1] * h[i] * h1*(h1 - 1.0)*(h1 - 1.0);
		}
		return spar;
	}

	void span(int spa)
	{
		string nnspan, zonename;
		double  temp, tem, cvl, qq2, qqw, qqr, a, t1;
		vector<double>  hx, hy, hz, hr, s, hu, hv, hw, hp, hpt, hh, hmiu, hwma;

		hx.resize(ny + 2);
		hy.resize(ny + 2);
		hz.resize(ny + 2);
		hr.resize(ny + 2);
		hu.resize(ny + 1);
		hw.resize(ny + 1);
		hp.resize(ny + 1);
		hpt.resize(ny + 1);
		hh.resize(ny + 1);
		hmiu.resize(ny + 1);
		hwma.resize(ny + 1);
		s.resize(ny + 2);
		hv.resize(ny + 2);

		vector<vector<double> >  hxx, hyy, hzz;

		hxx.resize(nz + 2);
		hyy.resize(nz + 2);
		hzz.resize(nz + 2);
		for (int i = 1; i<nz + 2; i++)
		{
			hxx[i].resize(nx + 2);
			hyy[i].resize(nx + 2);
			hzz[i].resize(nx + 2);
		}
		vector<vector<vector<double> > > huu, hvv, hww, hpp, hppt, hht, hwwma, hmmiu;

		huu.resize(nt + 1);
		hvv.resize(nt + 1);
		hww.resize(nt + 1);
		hpp.resize(nt + 1);
		hppt.resize(nt + 1);
		hht.resize(nt + 1);
		hwwma.resize(nt + 1);
		hmmiu.resize(nt + 1);
		for (int i = 1; i<nt + 1; i++)
		{
			huu[i].resize(nz + 1);
			hvv[i].resize(nz + 1);
			hww[i].resize(nz + 1);
			hpp[i].resize(nz + 1);
			hppt[i].resize(nz + 1);
			hht[i].resize(nz + 1);
			hwwma[i].resize(nz + 1);
			hmmiu[i].resize(nz + 1);
			for (int j = 1; j<nz + 1; j++)
			{
				huu[i][j].resize(nx + 1);
				hvv[i][j].resize(nx + 1);
				hww[i][j].resize(nx + 1);
				hpp[i][j].resize(nx + 1);
				hppt[i][j].resize(nx + 1);
				hht[i][j].resize(nx + 1);
				hwwma[i][j].resize(nx + 1);
				hmmiu[i][j].resize(nx + 1);

			}
		}

		temp = double(spa) / 100.0;

		for (int i = 1; i<nz + 2; i++)
			for (int k = 1; k<nx + 2; k++)
			{
				s[1] = 0;
				for (int j = 1; j<ny + 2; j++)
				{
					hx[j] = x[i][j][k];
					hy[j] = y[i][j][k];
					hz[j] = z[i][j][k];
					hr[j] = sqrt(hy[j] * hy[j] + hz[j] * hz[j]);

					if (j>1)
					{
						s[j] = s[j - 1] + sqrt(pow(hx[j] - hx[j - 1], 2) + pow(hr[j] - hr[j - 1], 2));
					}
				}

				t1 = s[ny + 1] * temp;

				hxx[i][k] = wl(s, hx, t1, ny + 1, 1);
				hyy[i][k] = wl(s, hy, t1, ny + 1, 1);
				hzz[i][k] = wl(s, hz, t1, ny + 1, 1);
			}

		for (int i = 1; i<nz + 1; i++)
			for (int k = 1; k<nx + 1; k++)
			{
				s[1] = 0;

				for (int j = 1; j<ny + 1; j++)
				{

					hx[j] = xx0[i][j][k];
					hy[j] = yy0[i][j][k];
					hz[j] = zz0[i][j][k];
					hr[j] = sqrt(hy[j] * hy[j] + hz[j] * hz[j]);

					if (j>1)
					{
						s[j] = s[j - 1] + sqrt(pow(hx[j] - hx[j - 1], 2) + pow(hr[j] - hr[j - 1], 2));
					}
				}

				t1 = s[ny] * temp;

				for (n = 1; n<nt + 1; n++)
				{
					for (int j = 1; j<ny + 1; j++)
					{
						y1 = yy0[i][j][k];
						z1 = zz0[i][j][k];

						hu[j] = q12[n][i][j][k] / q11[n][i][j][k];
						hv[j] = q13[n][i][j][k] / q11[n][i][j][k];
						hw[j] = q14[n][i][j][k] / q11[n][i][j][k];

						wx = hu[j];
						wy = hv[j] + rpm*z1;
						wz = hw[j] - rpm*y1;

						qqw = wx*wx + wy*wy + wz*wz;
						qqr = rpm*rpm*(z1*z1 + y1*y1);
						qq2 = hu[j] * hu[j] + hv[j] * hv[j] + hw[j] * hw[j];

						hp[j] = 0.40*(q15[n][i][j][k] - 0.50*q11[n][i][j][k] * qq2);

						tem = hp[j] / (q11[n][i][j][k] * rg);
						cvl = cvl0*pow((tem / t0), 1.5)*(t0 + ts) / (tem + ts);

						a = 1.40*hp[j] / q11[n][i][j][k];

						hwma[j] = sqrt(qqw / a);
						hh[j] = (q15[n][i][j][k] + hp[j]) / q11[n][i][j][k];
						hpt[j] = hp[j] * pow(1.0 - 0.50*qq2 / hh[j], -3.5);
						hmiu[j] = q16[n][i][j][k] / cvl;
					}
					huu[n][i][k] = wl(s, hu, t1, ny, 1);
					hvv[n][i][k] = wl(s, hv, t1, ny, 1);
					hww[n][i][k] = wl(s, hw, t1, ny, 1);
					hpp[n][i][k] = wl(s, hp, t1, ny, 1);
					hppt[n][i][k] = wl(s, hpt, t1, ny, 1);
					hht[n][i][k] = wl(s, hh, t1, ny, 1);
					hwwma[n][i][k] = wl(s, hwma, t1, ny, 1);
					hmmiu[n][i][k] = wl(s, hmiu, t1, ny, 1);
				}
			}


		//	*****************输出

		ostringstream os;
		os << spa;
		nnspan = os.str();


		fstream f21;
		string file = trim(nnspan) + "%span-" + trim(id_m) + "myid.dat";

		f21.open(file.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
		f21 << "VARIABLES=x,y,z,u,v,w,wma,pressure,pt,zht,ut";

		for (int n = 1; n < nt + 1; n++)
		{
			os.str("");
			os << n;
			zonename = os.str();
			zonename = "out_" + trim(zonename);

			int count = 0;

			if (n == 1)
			{
				f21 << "ZONE T=" << zonename << "I=" << nx + 1 << "J=" << ny + 1 << "K=" << nz + 1
					<< "DATAPACKING=BLOCK, VARLOCATION=([4-10]=CELLCENTERED),SOLUTIONTIME="
					<< setw(8) << setprecision(4) << double(n - 1) / double(nt);

				for (int i = 1; i<nz + 2; i++)
					for (int k = 1; k<nx + 2; k++)
					{
						f21 << setprecision(15) << hxx[i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
						//					if(myid==0)cout<<i<<' '<<k<<endl;
					}

				count = 0;

				for (int i = 1; i<nz + 2; i++)
					for (int k = 1; k<nx + 2; k++)
					{
						f21 << setprecision(15) << hyy[i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
					}

				count = 0;

				for (int i = 1; i<nz + 2; i++)
					for (int k = 1; k<nx + 2; k++)
					{
						f21 << setprecision(15) << hzz[i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
					}

				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int k = 1; k<nx + 1; k++)
					{
						f21 << setprecision(15) << huu[n][i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
					}

				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int k = 1; k<nx + 1; k++)
					{
						f21 << setprecision(15) << hvv[n][i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
					}

				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int k = 1; k<nx + 1; k++)
					{
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
						f21 << setprecision(15) << hww[n][i][k] << "      ";
						count = 0;
					}

				for (int i = 1; i<nz + 1; i++)
					for (int k = 1; k<nx + 1; k++)
					{
						f21 << setprecision(15) << hwwma[n][i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
					}

				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int k = 1; k<nx + 1; k++)
					{
						f21 << setprecision(15) << hpp[n][i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
					}

				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int k = 1; k<nx + 1; k++)
					{
						f21 << setprecision(15) << hppt[n][i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
					}

				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int k = 1; k<nx + 1; k++)
					{
						f21 << setprecision(15) << hht[n][i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
					}

				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int k = 1; k<nx + 1; k++)
					{
						f21 << setprecision(15) << hmmiu[n][i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
					}
			}
			else
			{
				f21 << "ZONE T=" << zonename << "I=" << nx + 1 << "J=" << ny + 1 << "K=" << nz + 1
					<< "DATAPACKING=BLOCK, VARSHARELIST=([1-3]=1),VARLOCATION=([4-10]=CELLCENTERED),SOLUTIONTIME="
					<< setw(8) << setprecision(4) << double(n - 1) / double(nt);

				for (int i = 1; i<nz + 1; i++)
					for (int k = 1; k<nx + 1; k++)
					{
						f21 << setprecision(15) << huu[n][i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
					}

				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int k = 1; k<nx + 1; k++)
					{
						f21 << setprecision(15) << hvv[n][i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
					}

				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int k = 1; k<nx + 1; k++)
					{
						f21 << setprecision(15) << hww[n][i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
					}

				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int k = 1; k<nx + 1; k++)
					{
						f21 << setprecision(15) << hwwma[n][i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
					}

				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int k = 1; k<nx + 1; k++)
					{
						f21 << setprecision(15) << hpp[n][i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
					}

				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int k = 1; k<nx + 1; k++)
					{
						f21 << setprecision(15) << hppt[n][i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
					}

				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int k = 1; k<nx + 1; k++)
					{
						f21 << setprecision(15) << hht[n][i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
					}

				count = 0;

				for (int i = 1; i<nz + 1; i++)
					for (int k = 1; k<nx + 1; k++)
					{
						f21 << setprecision(15) << hmmiu[n][i][k] << "      ";
						count++;
						if (count == 3)
						{
							f21 << endl;
							count = 0;
						}
					}
			}
		}

		f21.close();
	}

	//输出

	void output()
	{
		int  spa;

		x = xf;
		y = yf;
		z = zf;

		inlout();
		lbout();

		spa = 30;
		span(spa);
		spa = 95;
		span(spa);
	}

	//FMG多重网格

	void fmg(int style)
	{
		double str,ed;
		for (nng = nng0; nng<ng + 1; nng++)
		{
			geo();        //计算网格参数
			ye();         //边界值
						  //Vcycle

			for (it00 = 1; it00<cg[nng] + 1; it00++)    //表示各阶段的迭代步数，其数值后来由文件中读入
			{
				if (myid == 0)
					cout << it00 << "  of   " << cg[nng] << "  fmg march1 start-------------------" << endl;

				nitt = nitt + 1;                   //总迭代次数，在ini2中初始化为0
				rmsm = -11;

				//Timelevel

				for (n = 1; n<nt + 1; n++)
				{
					geon(nng);
					yen(nng);

					for (int i = 1; i<nz + 1; i++)   //驱动源项P2h，多重网格粗网格迭代中
						for (int j = 1; j<ny + 1; j++)
							for (int k = 1; k<nx + 1; k++)
							{
								qp1[i][j][k] = 0;
								qp2[i][j][k] = 0;
								qp3[i][j][k] = 0;
								qp4[i][j][k] = 0;
								qp5[i][j][k] = 0;
							}

					tsd();
					//时间导数项
					march1();
					if (myid == 0) cout << "march1" << endl;                 //3步R-K推进
					residual();               //迭代中最大残差量的计算 ，判断收敛\
											  											  					 if(myid==0) cout<<"residual"<<endl;
					ppp();                    //求原始变量
					if (myid == 0) cout << "ppp" << endl;
					bc();                     //边界条件（周向数据传递，实现并行）
					if (myid == 0) cout << "bc" << endl;
					step();                   //当地时间步长
					if (myid == 0) cout << "step" << endl;
					if(nng==3)
					{
						str=MPI_Wtime();
						ddd();
						ed=MPI_Wtime();
						if(myid==0)
							cout<<"ddd_time:"<<ed-str<<endl;
					}                    //对流通量
					else
						ddd();
					if (myid == 0) cout << "ddd" << endl;
					if(nng==3)
					{
						str=MPI_Wtime();
						qqq();
						ed=MPI_Wtime();
						if(myid==0)
							cout<<"qqq_time:"<<ed-str<<endl;
					}                    //对流通量
					else
						qqq();
					if (myid == 0) cout << "qqq" << endl;
					qqqv();                  //粘性通量
					if (myid == 0) cout << "qqqv" << endl;

					for (int i = 1; i<nz + 1; i++)
						for (int j = 1; j<ny + 1; j++)
							for (int k = 1; k<nx + 1; k++)
							{
								q01[i][j][k] = q11[n][i][j][k];
								q02[i][j][k] = q12[n][i][j][k];
								q03[i][j][k] = q13[n][i][j][k];
								q04[i][j][k] = q14[n][i][j][k];
								q05[i][j][k] = q15[n][i][j][k];
							}

					timl = cfl*0.125;
					pred(0);                //每步R-K推进方法

					if (nng == ng)
					{
						marchsa();            //湍流模型计算
						if (myid == 0) cout << "marchsa" << endl;
					}

					for (int n = 1; n<nt + 1; n++)
						for (int k = 1; k<nz + 1; k++)
							for (int j = 1; j<ny + 1; j++)
								for (int i = 1; i<nx + 1; i++)
								{
									q31[n][k][j][i] = q11[n][k][j][i];
									q32[n][k][j][i] = q12[n][k][j][i];
									q33[n][k][j][i] = q13[n][k][j][i];
									q34[n][k][j][i] = q14[n][k][j][i];
									q35[n][k][j][i] = q15[n][k][j][i];
									q36[n][k][j][i] = q16[n][k][j][i];
								}

					if (style == 1)
					{
						//********************粗网格计算*********************
						time_t = MPI_Wtime();
						for (ign = nng - 1; ign >= 1; ign--)
						{
							geon(ign);              //限制算子
							init();
							yen(ign);
							copr();                   //粗网格驱动项
							march2();                 //粗网格3步R-K推进
							update();                  //插值
							ppp();
							bc();
							step();
							dddc();
							qqq();
							qqqv();

							for (int i = 1; i<nz + 1; i++)
								for (int j = 1; j<ny + 1; j++)
									for (int k = 1; k<nx + 1; k++)
									{
										q01[i][j][k] = q11[n][i][j][k];
										q02[i][j][k] = q12[n][i][j][k];
										q03[i][j][k] = q13[n][i][j][k];
										q04[i][j][k] = q14[n][i][j][k];
										q05[i][j][k] = q15[n][i][j][k];
									}

							timl = cfl*0.1;
							pred(0);
							update();
						}

						//*****************粗网格计算完成*****************

						geon(nng);

						for (int i = 1; i<nt + 1; i++)
							for (int j = 1; j<nz + 1; j++)
								for (int k = 1; k<ny + 1; k++)
									for (int n = 1; n<nx + 1; n++)
									{
										q11[i][j][k][n] = q31[i][j][k][n];
										q12[i][j][k][n] = q32[i][j][k][n];
										q13[i][j][k][n] = q33[i][j][k][n];
										q14[i][j][k][n] = q34[i][j][k][n];
										q15[i][j][k][n] = q35[i][j][k][n];
										q16[i][j][k][n] = q36[i][j][k][n];
									}
					}
				}

				time_end = MPI_Wtime();

				f << setprecision(17) << nitt << ' ' << rmsm << ' ' << time_end - time_begin << endl;

				MPI_Allreduce(&rmsm, &rmsmmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

				if (rmsmmax == rmsm0)
					goto here;

				flow();     //出口流量

				if (myid == 0)
					cout << it00 << "  of   " << cg[nng] << "  fmg cycle ends-------------------" << endl;

			}                                         //end Vcycle
			if (nng<ng)
				update2();                              //插值跳阶段
			else
				output();                               //输出
		}
	here:;                                        //end FMGcycle
	}

}//cfd命名空间

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &cfd::myid);
	MPI_Comm_size(MPI_COMM_WORLD, &cfd::numprocs);

	cfd::time_begin = MPI_Wtime();

	cfd::ini();

	cfd::allocation();

	cfd::shu();

	cfd::ini2();

	cfd::dodd();

	cfd::distance();

	cfd::fmg(1);

	cfd::freealloc();

	MPI_Finalize();

	cout << "program finished" << endl;
	return 0;
}

