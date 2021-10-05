#include "stage35.h"
#include<stdio.h>
#include<math.h>
#include"mpi.h"
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include <sstream>
#include <string>
using namespace std;

namespace cfd {

	int m, n, mm, nn, mt, nt, nitt, ci, nci, nb;
	int i, j, k, l, ii, jj, kk, ifine, jfine, kfine;
	int nx, ny, nz, nx1, nx2, nxyp, ib, it, jb, jt, jt0;
	int ny0, iread, cc, xl, yl, zl, xll, xlln, xlnn, yln, zln;//�������������ʱ�����������������

	int *lb, *xln;

	double y1, z1, rr, sir, cor, vx, vy, vz, wx, wy, wz;
	double dim, en, pp, delt, cvl, a, qq2, two3, spa1, spa2;
	double c2, cfl, a4, a2, beta1, beta2, aa1, aa2, aa3;
	double val1, val2, val3, val4, val5, val6, val7, val8, val, dimin, ppin;

	double ta, timl, pt, ht, rout, pb0, pb1, period, rmsm;//��ѹ�����ʣ����ڣ���Ҷǰ��Ե��Ҷ����λ��
	double cvl0, t0, ts, cp, ccv, prt, prl, rg, cv1, cv2;
	double kap, sigmav, cb1, cb2, cw1, cw2, cw3, cr1, cr2, cr3;

	char nnspan[50], id_mm[50], id_sa[50];

	double *peb, *rpm, *ma, *temp0;

	double **hxx2, **hyy2, **hzz2;
	double **betax, **betay, **betaz;
	double **spa0, **spa01, **spa02, **hspa;
	double **hxx, **hyy, **hzz, **hxx1, **hyy1, **hzz1;

	double ***vv0;
	double ***dmini;
	double ***x, ***y, ***z;
	double ***x0, ***y0, ***z0;
	double ***xf0, ***yf0, ***zf0;
	double ***xx0, ***yy0, ***zz0;
	double ***xx01, ***yy01, ***zz01;
	double ***xx02, ***yy02, ***zz02;
	double ***xx03, ***yy03, ***zz03;
	double ***s1x, ***s1y, ***s1z;
	double ***s2x, ***s2y, ***s2z;
	double ***s3x, ***s3y, ***s3z;
	double ***sri, ***srj, ***srk;
	double ***pvx, ***pvy, ***pvz;
	double ***vth, ***vre, ***p, ***t, ***time;
	double ***q01, ***q02, ***q03, ***q04, ***q05, ***q06, ***gdf;// RKѭ��ǰ���غ����
	double ***q11, ***q12, ***q13, ***q14, ***q15, ***q16;//�غ����
	double ***ts1, ***ts2, ***ts3, ***ts4, ***ts5, ***ts6; // �����������ĺ�����
	double ***av1, ***av2, ***av3, ***av4, ***av5, ***av6; // �˹�ճ����
	double ***qc1, ***qc2, ***qc3, ***qc4, ***qc5, ***qc6; // ����ͨ��
	double ***qv1, ***qv2, ***qv3, ***qv4, ***qv5, ***qv6; // ճ��ͨ��

	double ****q21, ****q22, ****q23, ****q24, ****q25, ****q26; // ����ǰ��ʱ����������������õ�

//������
//char id_m[9], id_rm[9], id_lm[9], secname[9], id_num1[9];//���̺�,Ҷ�ź�,��Ҷ����ҶƬ��
	string id_m, id_rm, id_lm, secname, id_num1;//���̺�,Ҷ�ź�,��Ҷ����ҶƬ��
	double ****xxs, ****yys, ****zzs, **vv, ******v;
	int cli[3], bk[3], ek[3], bkk[3], ekk[3];
	int lm, lbm, slidm, rm, ssum;
	int myid, myidl, myidm, myidr, numprocs, numpp, ierr, rc;
	MPI_Status status;


#define pi 3.141592654

#define Imalloc1(array, num1)\
	array = (int*)malloc(sizeof(int) * num1)

#define Dmalloc1(array, num1)\
	array = (double*)malloc(sizeof(double) * num1)

	string& trim(string &s) {
		if (s.empty()) {
			return s;
		}
		s.erase(0, s.find_first_not_of(" "));
		s.erase(s.find_last_not_of(" ") + 1);
		return s;
	}

	void init(int argc, char **argv);
	void finished();
	int argv_process(int argc, char **argv);

	void ini();
	void allocation();
	void ini2();
	void distance();
	void geo();
	void spaa();
	void cfd();
	void slid();

	void overlapzone();
	void overlapgrid();
	void overlapadj();
	void march();
	void residual();
	void store();
	void probe();
	void test();
	void output0();
	void span(double spa11, double spa22);
	void chord();
	void wl(double *medx, double *medr, double *spax, int n1, int n2, double *spar);

	void init(int argc, char **argv) {
		MPI_Init(NULL, NULL);

		MPI_Comm_rank(MPI_COMM_WORLD, &myid);
		MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

		if(myid == 0) {
			printf("\n******************************\n");
			printf("         program begin\n");
			printf("******************************\n\n");
		}

		if(!argv_process(argc, argv)) {
			if(myid == 0) {
				printf("wrong argv\n");
			}
		}
	}

	void finished() {
		if(myid == 0) {
			printf("\n******************************\n");
			printf("       program finished\n");
			printf("******************************\n\n");
		}
		MPI_Finalize();
	}

	int argv_process(int argc, char **argv) {

		return 1;
	}

	/*
	 * ����˫ʱ�䲽�ƽ�������������ѭ��-����ʱ��ѭ��-����ʱ�����ѭ����
	 */
//int cfd() {
//	for(m = mm; m < mt + 1; m++) { //Periodcycle
//		if(m > mm)
//			nn = 1;
//		for(n = nn; n < nt + 1; n++) { //Physicaltimelevel �õ�������ʱ����ѭ��
//			//����Ҷ�����߲���ƽ�У��ʲ��ܼ򵥵�����ȴ�������ȣ�����ά����
//			if(slidm != 0) {
//				slid();//����Ҷ���ݣ���ͨ���õ�����Ҷ������ͨ���ı߽����񶥵�����ֵ
//				overlapzone();//�����ʵ���̱���Щ���ڽ��̲ü������ü���k���Ĳ���
//				overlapgrid();//��ǰ����ÿ�����񱻲�ֵ������
//				overlapadj();//�ص�����ռ��ǰ���̷�Χ������˵ ���ڽ���>���ص���Χ������bc�����������ķ�Χ�������ݽ������һ������x���
//			}
//			tsd();	//��������ֵ�ǰ��ʱ��ֵ
//			for(ci = 1; ci < nci + 1; ci++) { //Pseudotimelevel
//				nitt = nitt + 1;
//				march();		//����������1����1��(4��)
//				residual();		//�����������ϼ����ܶȲв����Ϊ����������ж�����
//				printf("%d, %lf, %d, %d, %d\n", nitt, rmsm, ci, n, m);//���ÿһ��������ĸ�ʱ�������в����¼ÿ����ǽ��ʱ��
//				if(rmsm > 3.0) {
//
//				} else if (rmsm <= 3.0) {
//
//				} else {
//					//�жϲв���ΪNAN����ֹͣ����
//					return 0;
//				}
//			}
//			store();			//�洢��ǰ��ǰһ���������ֵ���������м�ֵ
//			probe();
//			if(slidm != 0) {
//				test();
//			}
//		}
//		output0();
//	}
//	return 1;
//}

	/* ������ 1
	 * �����̶�����Ʋ���������
	 */
	void ini() {
		double temp, t1, t2;
		int nxm1, nym1, nzm1, num1, ml;
		int *nxm, *nym, *nzm, *ibm, *itm, *jbm, *jtm;
		cvl0 = 1.7161e-5;
		t0 = 273.16;	//288.15
		ts = 110.4;		//124
		rg = 287;
		cp = rg * 1.4 / 0.4;
		ccv = cp - rg;
		prl = 0.72;
		prt = 0.9;
		cv1 = 7.1;
		cv2 = 5;
		kap = 0.41;
		cb1 = 0.1355;
		cb2 = 0.622;
		sigmav = 2 / 3;
		cw1 = cb1 / (kap * kap) + (1 + cb2) / sigmav;
		cw2 = 0.3;
		cw3 = 2;
		cr1 = 1;
		cr2 = 2;	//09����Ϊ2��00����Ϊ12
		cr3 = 1;
		two3= 2 / 3;

		ifstream file10("ini3ji.dat");//��Ҷ������ʣ���ѹ���м伶ѹ��pb0�����ھ�ѹpb1
		//c2,�����ճ�����ӡ�mt�������ڸ���,nci��ĳһʱ��������ʱ������
		file10 >> iread >> nt >> beta1 >> beta2 >> cfl >> a2 >> a4 >> ht >> pt >> pb1 >> c2 >> nb >> mt >> nci >> spa1 >> spa2;

		Imalloc1(nxm, nb + 1);
		Imalloc1(nym, nb + 1);
		Imalloc1(nzm, nb + 1);
		Imalloc1(ibm, nb + 1);
		Imalloc1(itm, nb + 1);
		Imalloc1(jbm, nb + 1);
		Imalloc1(jtm, nb + 1);
		Imalloc1(lb,  nb + 1);
		Dmalloc1(rpm, nb + 1);
		Dmalloc1(ma,  nb + 1);

		for(l = 1; l < nb + 1; l++) {
			string str_temp;
			getline(file10, str_temp);
			getline(file10, str_temp);//����,д���в�����Ч
			//itm,jtm����̱ں�ĵ�һ��������������Ҫ��1
			file10 >> nxm[l] >> nym[l] >> nzm[l] >> ibm[l] >> itm[l] >> jbm[l] >> jtm[l] >> lb[l] >> rpm[l] >> ma[l];
			if(myid == 0)
				cout<< nxm[l] << nym[l] << nzm[l] << ibm[l] << itm[l] << jbm[l] << jtm[l] << lb[l] << rpm[l] << ma[l] <<endl;
			rpm[l] = rpm[l] * pi / 30;// ��λ�����rad/s--��ת���ٶȣ���������
		}

		file10.close();
		period = 2 * pi / abs(rpm[1]) / lb[1];   //����ҶҶƬ����Ϊ��ͨ����Ŀǰ�о����ĸ�Ҷ�ŵ�ҶƬ������ࣩ
		delt = period / double(nt);				//��ʵ����ʱ�䲽��
		numpp = 10;

		if(myid < numpp * lb[1])
			rm = 1;    //����������Ҷ�ź�1,2
		else
			rm = 2;

		ostringstream os1, os2;
		os1 << myid;
		id_m = os1.str();
		os2 << rm;
		id_rm = os2.str();

		nx = nxm[rm];
		ny = nym[rm];
		nz = nzm[rm];
		ib = ibm[rm];
		it = itm[rm] - 1;
		jb = jbm[rm];
		jt = jtm[rm] - 1;
		jt0 = jt;

		xf0 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		yf0 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		zf0 = malloc3D<double>(nx + 2, ny + 2, nz + 2);

		string file55_name = "grid-" + trim(id_rm) + ".dat";
		if(myid == 0)
			cout<< file55_name << endl;
		ifstream file55(file55_name.c_str());
		file55 >> nxm1 >> nym1 >> nzm1;

		//�����̶�����������
		for (k = 1; k < nzm1 + 1; k++)
			for (j = 1; j < nym1 + 1; j++)
				for (i = 1; i < nxm1 + 1; i++)
					file55 >> xf0[k][j][i] >> yf0[k][j][i] >> zf0[k][j][i];

		if(myid == 0)
			cout<< xf0[nzm1][nym1][nxm1];

		file55.close();

		/*��Ҷ���ֿ�*/
		xlnn = 3;      //x�����������������
		Imalloc1(xln, xlnn);

		if(rm == 1) {
			xln[0] = 2;//x����ÿ�������ڲ�����
			xln[2] = 1;
		} else {
			xln[0] = 1;
			xln[2] = 2;
		}

		xln[1] = 2;
		yln = 2;
		zln = 1;

		num1 = myid % numpp;   //ÿ��Ҷ���ڵĽ��̱��0~79

		if(num1 < xln[0] * yln * zln)
			xlln = num1 / (xln[0] * yln * zln);//x���������Ĵ��Ա��0��1��2
		else
			xlln = (num1 - xln[0] * yln * zln) / (xln[1] * yln * zln) + 1;

		ml = num1 % (xln[xlln] * zln);  //��ǰ������ÿ��y�ϵ�8�����̣�0~7����������4�����̣�0~3
		xl = ml / zln;   //�������ڲ����0~3��0~3��0~1
		xll = xl;//�ۺ���������£�������x�����ȫ�ֱ�ţ�0~9

		for(i = 1; i < xlln + 1; i++) {
			xll = xll + xln[i - 1];
		}

		for(i = 1; i < xlln + 1; i++) {
			num1 = num1 - xln[i - 1] * yln * zln;
		}

		num1 = num1 % (xln[xlln] * yln * zln);//ÿ������ڵĽ��̱��0~31��0~31��0~15
		yl = num1 / (xln[xlln] * zln);//y������
		zl = myid % zln;//z������

		/*����������ֵ*/
		nx = (it - ib + 1) / xln[1];
		/*��ǰ���̵�������*/
		if(rm == 1) {
			nx1 = nx;//��ǰ����ǰ��ÿ����������
			if(xll == xln[0] + xln[1] + xln[2] - 1)
				nx = 24;
		} else {
			nx1 = 24;//��ǰ����ǰ��ÿ����������
			if(xll == 0)
				nx = 24;
		}

		ny = ny / yln;
		nz = nz / zln;
		x0 = malloc3D<double>(nx + 3, ny + 2, nz + 2);
		y0 = malloc3D<double>(nx + 3, ny + 2, nz + 2);
		z0 = malloc3D<double>(nx + 3, ny + 2, nz + 2);

		if(rm == 1) {
			for(int k = 1; k < nz + 2; k++)
				for(int j = 1; j < ny + 2; j++)
					for(int i = 1; i < nx + 2; i++) {
						x0[k][j][i] = xf0[zl * nz + k][yl * ny + j][xll * nx1 + i];
						y0[k][j][i] = yf0[zl * nz + k][yl * ny + j][xll * nx1 + i];
						z0[k][j][i] = zf0[zl * nz + k][yl * ny + j][xll * nx1 + i];
					}
		} else {
			for(int k = 1; k < nz + 2; k++)
				for(int j = 1; j < ny + 2; j++)
					for(int i = 1; i < nx + 2; i++) {
						x0[k][j][i] = xf0[zl * nz + k][yl * ny + j][xll * nx - nx + nx1 + i];
						y0[k][j][i] = yf0[zl * nz + k][yl * ny + j][xll * nx - nx + nx1 + i];
						z0[k][j][i] = zf0[zl * nz + k][yl * ny + j][xll * nx - nx + nx1 + i];
					}
		}

		//��ת���Ƶ���Ȧ����
		x = malloc3D<double>(nx + 3, ny + 2, nz + 2);
		y = malloc3D<double>(nx + 3, ny + 2, nz + 2);
		z = malloc3D<double>(nx + 3, ny + 2, nz + 2);

		lbm = myid / numpp;//���������ڵ�ҶƬ��0,1,,35,36,,81
		if(lbm < lb[1])
			lm = lbm;
		else
			lm = lbm - lb[1];//������󣬸���������Ҷ�ŵ�ҶƬ�ţ�0,,35.0,,45

		temp = 2 / double(lb[rm]) * pi * lm;
		cor = cos(temp);								//0             ^
		sir = sin(temp);								//1             ^        ҶƬ������ת
		for(k = 1; k < nz + 2; k++) { 					//nz            ^
			for(j = 1; j < ny + 2; j++) { 				//nz+1          1
				for(i = 1; i < nx + 2; i++) {
					x[k][j][i] = x0[k][j][i];
					y[k][j][i] = y0[k][j][i] * cor + z0[k][j][i] * sir;//��ҶƬ��ת������ת,���෴����Ϊ����0��nz+1���꣬1-nz��˳����ҶƬ��ת�����෴
					z[k][j][i] = - y0[k][j][i] * sir + z0[k][j][i] * cor;
				}
			}
		}

		/*ҶƬ��y����Χ*/
		if(yln == 2) {
			if(yl == 0)
				jt = ny;
			else
				jt = jt - ny;
		}

		/*�ж��Ƿ�ִ�л�������*/
		slidm = 0;
		if(rm == 1) {
			if(xll == xln[0] + xln[1] + xln[2] - 1)
				slidm = 1;
		} else {
			if(xll == 0)
				slidm = 2;
		}
	}

///* ��v���ȷ��
// * �����ڴ�
// */
	void allocation() {
		gdf = malloc3D<double>(16, ny + 2, nz + 2);

		q11 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		q12 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		q13 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		q14 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		q15 = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		q16 = malloc3D<double>(nx + 2, ny + 2, nz + 2);

		q21 = malloc4D<double>(nx + 1, ny + 1, nz + 1, 3);
		q22 = malloc4D<double>(nx + 1, ny + 1, nz + 1, 3);
		q23 = malloc4D<double>(nx + 1, ny + 1, nz + 1, 3);
		q24 = malloc4D<double>(nx + 1, ny + 1, nz + 1, 3);
		q25 = malloc4D<double>(nx + 1, ny + 1, nz + 1, 3);
		q26 = malloc4D<double>(nx + 1, ny + 1, nz + 1, 3);

		dmini = malloc3D<double>(nx + 1, ny + 1, nz + 1);

		hxx = malloc2D<double>(nx + 2, nz + 2);
		hyy = malloc2D<double>(nx + 2, nz + 2);
		hzz = malloc2D<double>(nx + 2, nz + 2);

		spa0  = malloc2D<double>(nx + 1, nz + 1);

		hxx1 = malloc2D<double>(nx + 2, nz + 2);
		hyy1 = malloc2D<double>(nx + 2, nz + 2);
		hzz1 = malloc2D<double>(nx + 2, nz + 2);

		hxx2 = malloc2D<double>(nx + 2, nz + 2);
		hyy2 = malloc2D<double>(nx + 2, nz + 2);
		hzz2 = malloc2D<double>(nx + 2, nz + 2);

		spa01 = malloc2D<double>(nx + 1, nz + 1);
		spa02 = malloc2D<double>(nx + 1, nz + 1);

		hspa  = malloc2D<double>(nx + 1, nz + 1);

		xxs = malloc4D<double>(4, ny + 2, nz + 2, lb[1] + lb[2]);//���յ�������Ҷ�ű߽�����
		yys = malloc4D<double>(4, ny + 2, nz + 2, lb[1] + lb[2]);
		zzs = malloc4D<double>(4, ny + 2, nz + 2, lb[1] + lb[2]);

		vv = malloc2D<double>(ny + 1, nz + 1);

		//allocate(v, ny + 1, nz + 1, 3, ny + 1, nz + 1, 4);//ĳ���������ص���ռ��ǰ����������

		betax = malloc2D<double>(ny + 1, nz + 1);
		betay = malloc2D<double>(ny + 1, nz + 1);
		betaz = malloc2D<double>(ny + 1, nz + 1);

		Dmalloc1(peb, ny + 1);//���ھ�ѹ

		ts1 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		ts2 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		ts3 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		ts4 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		ts5 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		ts6 = malloc3D<double>(nx + 1, ny + 1, nz + 1);

		s1x = malloc3D<double>(nx + 1, ny + 1, nz + 3);//z����
		s1y = malloc3D<double>(nx + 1, ny + 1, nz + 3);
		s1z = malloc3D<double>(nx + 1, ny + 1, nz + 3);
		s2x = malloc3D<double>(nx + 3, ny + 1, nz + 1);//x����
		s2y = malloc3D<double>(nx + 3, ny + 1, nz + 1);
		s2z = malloc3D<double>(nx + 3, ny + 1, nz + 1);
		s3x = malloc3D<double>(nx + 1, ny + 3, nz + 1);//y����
		s3y = malloc3D<double>(nx + 1, ny + 3, nz + 1);
		s3z = malloc3D<double>(nx + 1, ny + 3, nz + 1);

		vv0 = malloc3D<double>(nx + 1, ny + 1, nz + 1);

		xx01 = malloc3D<double>(nx + 1, ny + 1, nz + 2);
		yy01 = malloc3D<double>(nx + 1, ny + 1, nz + 2);
		zz01 = malloc3D<double>(nx + 1, ny + 1, nz + 2);
		xx02 = malloc3D<double>(nx + 2, ny + 1, nz + 1);
		yy02 = malloc3D<double>(nx + 2, ny + 1, nz + 1);
		zz02 = malloc3D<double>(nx + 2, ny + 1, nz + 1);
		xx03 = malloc3D<double>(nx + 1, ny + 2, nz + 1);
		yy03 = malloc3D<double>(nx + 1, ny + 2, nz + 1);
		zz03 = malloc3D<double>(nx + 1, ny + 2, nz + 1);

		xx0 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		yy0 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		zz0 = malloc3D<double>(nx + 1, ny + 1, nz + 1);

		q01 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		q02 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		q03 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		q04 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		q05 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		q06 = malloc3D<double>(nx + 1, ny + 1, nz + 1);

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

		time = malloc3D<double>(nx + 1, ny + 1, nz + 1);

		sri = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		srj = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		srk = malloc3D<double>(nx + 2, ny + 2, nz + 2);
	}


	fstream file4, file7, file8, file11, file20, file121, file122, file123, file124, file125, file126, file131, file132;
	fstream file133, file134, file231, file232, file233, file234, file300;
	/*
	 * ��ʼ������
	 */
	void ini2() {
		double vxx, vrr, vtt, vee, vth1, vre1;
		int num1, num2, num3;
		/* ȷ������б�ʽ� */
		vxx = 1;
		vtt = vxx * beta1;//չ��������*vxx
		vrr = vxx * beta2;//����������*vxx
		vee = sqrt(vxx * vxx + vrr * vrr + vtt * vtt);
		for (k = 1; k < nz + 1; k++)
			for (j = 1; j < ny + 1; j++) {
				y1 = 0.25 * (y[k][j][1] + y[k][j + 1][1] + y[k + 1][j][1] + y[k + 1][j + 1][1]);
				z1 = 0.25 * (z[k][j][1] + z[k][j + 1][1] + z[k + 1][j][1] + z[k + 1][j + 1][1]);

				rr = sqrt(y1 * y1 + z1 * z1);

				sir = z1 / rr;
				cor = y1 / rr;

				for (int i = 1; i < nt + 1; i++) {
					betax[k][j] = vxx / vee;//!��һ����Ҷ
					betay[k][j] = (vrr * cor - vtt * sir) / vee;
					betaz[k][j] = (vrr * sir + vtt * cor) / vee;
				}
			}

		rout = 3.5 * pt / ht;
		dim = rout * pow(1 + 0.2 * ma[rm] * ma[rm], -2.5);
		pp = pt * pow(1 + 0.2 * ma[rm] * ma[rm], -3.5);
		a = sqrt(1.4 * pp / dim);
		en = 2.5 * pp + 0.5 * dim * ma[rm] * ma[rm] * a * a;//�ļ������Ǿ��������
		dimin = dim;
		ppin = pp;

		for (k = 1; k < nz + 1; k++)
			for (j = 1; j < ny + 1; j++)
				for (i = 1; i < nx + 1; i++) {
					q11[k][j][i] = dim;
					q12[k][j][i] = dim * a * ma[rm];
					q15[k][j][i] = en;
					q16[k][j][i] = 200 * cvl0;
				}


		if(iread == 0) {
			num1 = myid % numpp;
			if(rm == 1) {

			} else {
				num1 = num1 + 10;
			}
			ostringstream os;
			os << num1;
			id_num1 = os.str();

			string file300_name = "conservative_var-" + trim(id_num1) + "myid.dat";//ʹ�ö������������һ��������ѹ������Ϊ����
			ifstream file300(file300_name.c_str());
			//�˴�����Ҫ�ص��ļ���ʼ
			//file300.clear();
			//file300.seekg(ios::beg);

			file300 >> num2 >> num3;//������������Ľ׶ε������

			for (k = 1; k < nz + 1; k++)
				for (j = 1; j < ny + 1; j++)
					for (i = 1; i < nx + 1; i++)
						file300 >> q11[k][j][i];

			for (k = 1; k < nz + 1; k++)
				for (j = 1; j < ny + 1; j++)
					for (i = 1; i < nx + 1; i++)
						file300 >> q12[k][j][i];

			for (k = 1; k < nz + 1; k++)
				for (j = 1; j < ny + 1; j++)
					for (i = 1; i < nx + 1; i++)
						file300 >> q13[k][j][i];

			for (k = 1; k < nz + 1; k++)
				for (j = 1; j < ny + 1; j++)
					for (i = 1; i < nx + 1; i++)
						file300 >> q14[k][j][i];

			for (k = 1; k < nz + 1; k++)
				for (j = 1; j < ny + 1; j++)
					for (i = 1; i < nx + 1; i++)
						file300 >> q15[k][j][i];

			for (k = 1; k < nz + 1; k++)
				for (j = 1; j < ny + 1; j++)
					for (i = 1; i < nx + 1; i++)
						file300 >> q16[k][j][i];

			file300.close();

			for(k = 1; k < nz + 1; k++) { //0Ҷ���ϵ�q13,q14ת��Ϊ��ǰҶ��
				for(j = 1; j < ny + 1; j++) {
					for(i = 1; i < nx + 1; i++) {
						vy = q13[k][j][i] / q11[k][j][i];
						vz = q14[k][j][i] / q11[k][j][i];
						y1 = 0.125 * (y0[k][j][i] + y0[k + 1][j][i] + y0[k][j][i + 1] + y0[k + 1][j][i + 1]
						              + y0[k][j + 1][i] + y0[k + 1][j + 1][i] + y0[k][j + 1][i + 1] + y0[k + 1][j + 1][i + 1]);
						z1 = 0.125 * (z0[k][j][i] + z0[k + 1][j][i] + z0[k][j][i + 1] + z0[k + 1][j][i + 1]
						              + z0[k][j + 1][i] + z0[k + 1][j + 1][i] + z0[k][j + 1][i + 1] + z0[k + 1][j + 1][i + 1]);
						rr = sqrt(y1 * y1 + z1 * z1);
						sir = z1 / rr;
						cor = y1 / rr;
						vth1 = vz * cor - vy * sir;//Ҷ��y�ٶ�  ת��Ϊ0Ҷ���ϵ�Բ���ٶ�
						vre1 = vz * sir + vy * cor;//����z�ٶ�

						y1 = 0.125 * (y[k][j][i] + y[k + 1][j][i] + y[k][j][i + 1] + y[k + 1][j][i + 1]
						              + y[k][j + 1][i] + y[k + 1][j + 1][i] + y[k][j + 1][i + 1] + y[k + 1][j + 1][i + 1]);//Ҷ��y�ٶ�  ת��Ϊ��ǰҶ���ϵ�Բ���ٶ�
						z1 = 0.125 * (z[k][j][i] + z[k + 1][j][i] + z[k][j][i + 1] + z[k + 1][j][i + 1]
						              + z[k][j + 1][i] + z[k + 1][j + 1][i] + z[k][j + 1][i + 1] + z[k + 1][j + 1][i + 1]);
						rr = sqrt(y1 * y1 + z1 * z1);
						sir = z1 / rr;
						cor = y1 / rr;
						vy = vre1 * cor - vth1 * sir;
						vz = vre1 * sir + vth1 * cor;
						q13[k][j][i] = q11[k][j][i] * vy;
						q14[k][j][i] = q11[k][j][i] * vz;
					}
				}
			}

			mm = 1;
			nn = 2;
			nitt = 0;
			for(i = 1; i <= 2; i++) {
				for(int kter = 1; kter < nz + 1; kter++) {
					for(int jter = 1; jter < ny + 1; jter++) {
						for(int iter = 1; iter < nx + 1; iter++) {
							q21[i][kter][jter][iter] = q11[kter][jter][iter];//�м�����������ʱ��Ҫ֪��ǰ����ʱ��ֵ
							q22[i][kter][jter][iter] = q12[kter][jter][iter];
							q23[i][kter][jter][iter] = q13[kter][jter][iter];
							q24[i][kter][jter][iter] = q14[kter][jter][iter];
							q25[i][kter][jter][iter] = q15[kter][jter][iter];
							q26[i][kter][jter][iter] = q16[kter][jter][iter];
						}
					}
				}
			}

			if(myid == 0) {
				file4.open("convergence-inflow.dat", std::fstream::out | std::fstream::trunc);
			}

			if(rm == 2 && myid == lb[1] * numpp + xln[0] * yln + xln[1] * yln + xln[2] - 1) {
				file7.open("convergence-outflow.dat", std::fstream::out | std::fstream::trunc);
				file8.open("pbi.dat", std::fstream::out | std::fstream::trunc);
				file11.open("eff.dat", std::fstream::out | std::fstream::trunc);
			}
			/*
			/*������ȳ�ʼ����error�ļ���д*/
			string file20_name = "error-" + trim(id_m) + ".dat";
			file20.open("convergence-outflow.dat", std::fstream::out | std::fstream::trunc);
			if(lm == 0 && slidm != 0 && yl == 0) {

				string file121_name = "flu" + trim(id_rm) + "-1.dat";
				string file122_name = "flu" + trim(id_rm) + "-2.dat";
				string file123_name = "flu" + trim(id_rm) + "-3.dat";
				string file124_name = "flu" + trim(id_rm) + "-4.dat";
				string file125_name = "flu" + trim(id_rm) + "-5.dat";
				string file126_name = "flu" + trim(id_rm) + "-6.dat";
				file121.open(file121_name.c_str(), std::fstream::out | std::fstream::trunc);
				file122.open(file122_name.c_str(), std::fstream::out | std::fstream::trunc);
				file123.open(file123_name.c_str(), std::fstream::out | std::fstream::trunc);
				file124.open(file124_name.c_str(), std::fstream::out | std::fstream::trunc);
				file125.open(file125_name.c_str(), std::fstream::out | std::fstream::trunc);
				file126.open(file126_name.c_str(), std::fstream::out | std::fstream::trunc);

			}
			if(rm == 1) { //��Ҷ����Ҷǰ��Եÿ��Ҷ������������̽��
				if(myid == 3 + lbm * numpp) {
					string file131_name = "rl-u-" + trim(id_m) + ".dat";
					string file132_name = "rl-p-" + trim(id_m) + ".dat";
					file131.open(file131_name.c_str(), std::fstream::out | std::fstream::trunc);
					file132.open(file132_name.c_str(), std::fstream::out | std::fstream::trunc);
				}
				if(myid == 9 + lbm * numpp) {
					string file133_name = "rt-u-" + trim(id_m) + ".dat";
					string file134_name = "rt-p-" + trim(id_m) + ".dat";
					file133.open(file133_name.c_str(), std::fstream::out | std::fstream::trunc);
					file134.open(file134_name.c_str(), std::fstream::out | std::fstream::trunc);
				}
			}
			if(rm == 2) {
				if(myid == 1 + lbm * numpp) {
					string file231_name = "sl-u-" + trim(id_m) + ".dat";
					string file232_name = "sl-p-" + trim(id_m) + ".dat";
					file231.open(file231_name.c_str(), std::fstream::out | std::fstream::trunc);
					file232.open(file232_name.c_str(), std::fstream::out | std::fstream::trunc);
				}
				if(myid == 8 + lbm * numpp) {
					string file233_name = "st-u-" + trim(id_m) + ".dat";
					string file234_name = "st-p-" + trim(id_m) + ".dat";
					file233.open(file233_name.c_str(), std::fstream::out | std::fstream::trunc);
					file234.open(file234_name.c_str(), std::fstream::out | std::fstream::trunc);
				}
			}
		} else if (iread == 1) {
//			//��Ϊ�������Ʋ�Ӱ���ȡ�����Ծ���300����
			string file300_name = "pause-" + trim(id_m) + "myid.dat";//unformattedʹ�ö����Ƹ�ʽ����
			ifstream file300(file300_name.c_str());
			//�˴�����Ҫ�ص��ļ���ʼ
			//file300.clear();
			//file300.seekg(ios::beg);

			file300 >> mm;            //�����������������ѭ����
			file300 >> nn;            //������������ĵ�n����ʱ��
			file300 >> nitt;          //������������ĵ�nitt������

			for(int iteration = 1; iteration <= 2; iteration++)
				for(k = 1; k < nz + 1; k++)
					for(j = 1; j < ny + 1; j++)
						for(i = 1; i < nx + 1; i++)
							file300 >> q21[iteration][k][j][i];

			for(int iteration = 1; iteration <= 2; iteration++)
				for(k = 1; k < nz + 1; k++)
					for(j = 1; j < ny + 1; j++)
						for(i = 1; i < nx + 1; i++)
							file300 >> q22[iteration][k][j][i];

			for(int iteration = 1; iteration <= 2; iteration++)
				for(k = 1; k < nz + 1; k++)
					for(j = 1; j < ny + 1; j++)
						for(i = 1; i < nx + 1; i++)
							file300 >> q23[iteration][k][j][i];

			for(int iteration = 1; iteration <= 2; iteration++)
				for(k = 1; k < nz + 1; k++)
					for(j = 1; j < ny + 1; j++)
						for(i = 1; i < nx + 1; i++)
							file300 >> q24[iteration][k][j][i];

			for(int iteration = 1; iteration <= 2; iteration++)
				for(k = 1; k < nz + 1; k++)
					for(j = 1; j < ny + 1; j++)
						for(i = 1; i < nx + 1; i++)
							file300 >> q25[iteration][k][j][i];

			for(int iteration = 1; iteration <= 2; iteration++)
				for(k = 1; k < nz + 1; k++)
					for(j = 1; j < ny + 1; j++)
						for(i = 1; i < nx + 1; i++)
							file300 >> q26[iteration][k][j][i];

			file300.close();
			for(k = 1; k < nz + 1; k++) {
				for(j = 1; j < ny + 1; j++) {
					for(i = 1; i < nx + 1; i++) {
						q11[k][j][i] = q21[2][k][j][i];
						q12[k][j][i] = q22[2][k][j][i];
						q13[k][j][i] = q23[2][k][j][i];
						q14[k][j][i] = q24[2][k][j][i];
						q15[k][j][i] = q25[2][k][j][i];
						q16[k][j][i] = q26[2][k][j][i];
					}
				}
			}
			nn = nn + 1;
			if(myid == 0) {
				file4.open("convergence-inflow.dat", std::fstream::out | std::fstream::app);
			}

			if(rm == 2 && myid == lb[1] * numpp + xln[0] * yln + xln[1] * yln + xln[2] - 1) {
				file7.open("convergence-outflow.dat", std::fstream::out | std::fstream::app);
				file8.open("pbi.dat", std::fstream::out | std::fstream::app);
				file11.open("eff.dat", std::fstream::out | std::fstream::app);

			}
			file20.open("convergence-outflow.dat", std::fstream::out | std::fstream::app);
			if(lm == 0 && slidm != 0 && yl == 0) {
				string file121_name = "flu" + trim(id_rm) + "-1.dat";
				string file122_name = "flu" + trim(id_rm) + "-2.dat";
				string file123_name = "flu" + trim(id_rm) + "-3.dat";
				string file124_name = "flu" + trim(id_rm) + "-4.dat";
				string file125_name = "flu" + trim(id_rm) + "-5.dat";
				string file126_name = "flu" + trim(id_rm) + "-6.dat";
				file121.open(file121_name.c_str(), std::fstream::out | std::fstream::app);
				file122.open(file122_name.c_str(), std::fstream::out | std::fstream::app);
				file123.open(file123_name.c_str(), std::fstream::out | std::fstream::app);
				file124.open(file124_name.c_str(), std::fstream::out | std::fstream::app);
				file125.open(file125_name.c_str(), std::fstream::out | std::fstream::app);
				file126.open(file126_name.c_str(), std::fstream::out | std::fstream::app);
			}

			if(rm == 1) { //��Ҷ����Ҷǰ��Եÿ��Ҷ������������̽��
				if(myid == 3 + lbm * numpp) {
					string file131_name = "rl-u-" + trim(id_m) + ".dat";
					string file132_name = "rl-p-" + trim(id_m) + ".dat";
					file131.open(file131_name.c_str(), std::fstream::out | std::fstream::app);
					file132.open(file132_name.c_str(), std::fstream::out | std::fstream::app);
				}
				if(myid == 9 + lbm * numpp) {
					string file133_name = "rt-u-" + trim(id_m) + ".dat";
					string file134_name = "rt-p-" + trim(id_m) + ".dat";
					file133.open(file133_name.c_str(), std::fstream::out | std::fstream::app);
					file134.open(file134_name.c_str(), std::fstream::out | std::fstream::app);
				}
			}
			if(rm == 2) {
				if(myid == 1 + lbm * numpp) {
					string file231_name = "sl-u-" + trim(id_m) + ".dat";
					string file232_name = "sl-p-" + trim(id_m) + ".dat";
					file231.open(file231_name.c_str(), std::fstream::out | std::fstream::app);
					file232.open(file232_name.c_str(), std::fstream::out | std::fstream::app);
				}
				if(myid == 8 + lbm * numpp) {
					string file233_name = "st-u-" + trim(id_m) + ".dat";
					string file234_name = "st-p-" + trim(id_m) + ".dat";
					file233.open(file233_name.c_str(), std::fstream::out | std::fstream::app);
					file234.open(file234_name.c_str(), std::fstream::out | std::fstream::app);
				}
			}
		}
	}

	void distance() {
		int iil, iir, kkl, kkr, iib, iit, jjb, jjt, nx3, nx4, jtt, ny3, ny4;
		double d, dy1, dy2, dy, dz1, dz2, dz;
		double **xwu, **ywu, **zwu, **xwd, **ywd, **zwd, **xwf, **ywf, **zwf, **xwb;
		double **ywb, **zwb; //����Ҷ�ŵı�������
		double **xd, **yd, **zd, **xu, **yu, **zu, **xb, **yb, **zb, **xf, **yf, **zf;
		double ***xx00, ***yy00, ***zz00;


		/*****���������ĵ�����*****/
		xx00 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		yy00 = malloc3D<double>(nx + 1, ny + 1, nz + 1);
		zz00 = malloc3D<double>(nx + 1, ny + 1, nz + 1);

		for(k = 1; k < nz + 1; k++) {
			for(int j = 1; j < ny + 1; j++) {
				for(int i = 1; i < nx + 1; i++) {
					xx00[k][j][i] = 0.125 * (x[k][j][i] + x[k + 1][j][i] + x[k][j][i + 1] + x[k + 1][j][i + 1]
					                         + x[k][j + 1][i] + x[k + 1][j + 1][i] + x[k][j + 1][i + 1] + x[k + 1][j + 1][i + 1]);
					yy00[k][j][i] = 0.125 * (y[k][j][i] + y[k + 1][j][i] + y[k][j][i + 1] + y[k + 1][j][i + 1]
					                         + y[k][j + 1][i] + y[k + 1][j + 1][i] + y[k][j + 1][i + 1] + y[k + 1][j + 1][i + 1]);
					zz00[k][j][i] = 0.125 * (z[k][j][i] + z[k + 1][j][i] + z[k][j][i + 1] + z[k + 1][j][i + 1]
					                         + z[k][j + 1][i] + z[k + 1][j + 1][i] + z[k][j + 1][i + 1] + z[k + 1][j + 1][i + 1]);
				}
			}
		}

		nxyp = (it - ib + 1) / xln[1];

		/*��������*/
		/*yҶ�߷������±�������������*/
		xwd = malloc2D<double>(nx + 1, nz + 1);//j=1
		ywd = malloc2D<double>(nx + 1, nz + 1);
		zwd = malloc2D<double>(nx + 1, nz + 1);
		xwu = malloc2D<double>(nx + 1, nz + 1);//j=ny
		ywu = malloc2D<double>(nx + 1, nz + 1);
		zwu = malloc2D<double>(nx + 1, nz + 1);

		//�������븺���±������(xd, yd, zd, xu, yu, zu),�˴���temp1����Ϊ��̬���� ,ԭ��ΧΪ(-nxyp:2*nxyp,1:nz)
		double **xd_temp1, **yd_temp1, **zd_temp1, **xu_temp1, **yu_temp1, **zu_temp1;
		double *xd_temp2[nz + 1], *yd_temp2[nz + 1], *zd_temp2[nz + 1], *xu_temp2[nz + 1], *yu_temp2[nz + 1], *zu_temp2[nz + 1];

		//��������ռ�Ķ�ά���飬ע����Ǵ˴�ʹ�õķ�ʽ���ɵ������෴�ģ����ɾ�̬Ϊ xd_temp1[nz + 1][3 * nxyp + 1]
		xd_temp1 = malloc2D<double>(3 * nxyp + 1, nz + 1);//j=1
		yd_temp1 = malloc2D<double>(3 * nxyp + 1, nz + 1);
		zd_temp1 = malloc2D<double>(3 * nxyp + 1, nz + 1);
		xu_temp1 = malloc2D<double>(3 * nxyp + 1, nz + 1);//j=ny
		yu_temp1 = malloc2D<double>(3 * nxyp + 1, nz + 1);
		zu_temp1 = malloc2D<double>(3 * nxyp + 1, nz + 1);
//		//�����ά��
		for(i = 0; i < nz + 1; i++) {
			xd_temp2[i] = xd_temp1[i] + nxyp;
			yd_temp2[i] = yd_temp1[i] + nxyp;
			zd_temp2[i] = zd_temp1[i] + nxyp;
			xu_temp2[i] = xu_temp1[i] + nxyp;
			yu_temp2[i] = yu_temp1[i] + nxyp;
			zu_temp2[i] = zu_temp1[i] + nxyp;
		}
//		//�����ά��
		xd = xd_temp2;
		yd = yd_temp2;
		zd = zd_temp2;
		xu = xu_temp2;
		yu = yu_temp2;
		zu = zu_temp2;

		if(myid == 0) {
			cout<<endl<<xln[xlln]<<endl;
		}
		if(yl == 0) {
			for(k = 1; k < nz + 1; k++)
				for(i = 1; i < nx + 1; i++) {
					xwd[k][i] = 0.25 * (x[k][1][i] + x[k][1][i + 1] + x[k + 1][1][i] + x[k + 1][1][i + 1]);
					ywd[k][i] = 0.25 * (y[k][1][i] + y[k][1][i + 1] + y[k + 1][1][i] + y[k + 1][1][i + 1]);
					zwd[k][i] = 0.25 * (z[k][1][i] + z[k][1][i + 1] + z[k + 1][1][i] + z[k + 1][1][i + 1]);
				}

			//�˴�������Ϊ������,���Բ�ȡ�ഫ���ݵķ�ʽ����
			//�˴����ôӵ�һ�п�ʼ���䣬����Ϊ������Ĵ��䷽ʽ����һ��

//				MPI_Sendrecv(&xwd[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 760,
//				             &xwu[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 860, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&ywd[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 770,
//				             &ywu[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 870, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&zwd[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 780,
//				             &zwu[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 880, MPI_COMM_WORLD, &status);
		} else if(yl == 1) {
			for(k = 1; k < nz + 1; k++) {
				for(i = 1; i < nx + 1; i++) {
					xwu[k][i] = 0.25 * (x[k][ny + 1][i] + x[k][ny + 1][i + 1] + x[k + 1][ny + 1][i] + x[k + 1][ny + 1][i + 1]);
					ywu[k][i] = 0.25 * (y[k][ny + 1][i] + y[k][ny + 1][i + 1] + y[k + 1][ny + 1][i] + y[k + 1][ny + 1][i + 1]);
					zwu[k][i] = 0.25 * (z[k][ny + 1][i] + z[k][ny + 1][i + 1] + z[k + 1][ny + 1][i] + z[k + 1][ny + 1][i + 1]);
				}
			}
//				MPI_Sendrecv(&xwu[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 760,
//				             &xwd[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 860, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&ywu[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 770,
//				             &ywd[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 870, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&zwu[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 780,
//				             &zwd[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 880, MPI_COMM_WORLD, &status);
		}

		if(xll > 0) {
			if(xl == 1) {
				j = myid - 1;
			} else {
				if(yl == 0)
					j = myid - xln[xlln - 1] - 1;
				else
					j = myid - xln[xlln] - 1;
			}
			if(rm == 2 && xll == 1)
				nx1 = 24;
			else
				nx1 = nxyp;

//				MPI_Sendrecv(&xwu[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 70,
//				             &xu[-nx1][0], nx1 * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 60, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&ywu[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 71,
//				             &yu[-nx1][0], nx1 * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 61, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&zwu[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 72,
//				             &zu[-nx1][0], nx1 * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 62, MPI_COMM_WORLD, &status);
//
//				MPI_Sendrecv(&xwd[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 73,
//				             &xd[-nx1][0], nx1 * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 63, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&ywd[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 74,
//				             &yd[-nx1][0], nx1 * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 64, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&zwd[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 75,
//				             &zd[-nx1][0], nx1 * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 65, MPI_COMM_WORLD, &status);
		} else {
			nx1 = - 1;
		}

		for (k = 1; k < nz + 1; k++) {
			for (i = 1; i < nx + 1; i++) {
				xd[k][i] = xwd[k][i];
				yd[k][i] = ywd[k][i];
				zd[k][i] = zwd[k][i];
				xu[k][i] = xwu[k][i];
				yu[k][i] = ywu[k][i];
				zu[k][i] = zwu[k][i];
			}
		}

		if(xll < xln[0] + xln[1] + xln[2] - 1) {
			if(xl == 0) {
				if(slidm == 2) {
					if(yl == 0)
						j = myid + xln[xlln] + 1;
					else
						j = myid + xln[xlln + 1] + 1;
				} else
					j = myid + 1;
			} else {
				if(yl == 0)
					j = myid + xln[xlln] + 1;
				else
					j = myid + xln[xlln + 1] + 1;
			}

			if(rm == 1 && xll == 3)
				nx2 = 24;
			else
				nx2 = nxyp;

//			if(myid + xln[xlln] < 820) {//��ȷ���Բ���
//				MPI_Sendrecv(&xwu[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 60,
//				             &xu[nx + 1][0], nx2 * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 70, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&ywu[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 61,
//				             &yu[nx + 1][0], nx2 * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 71, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&zwu[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 62,
//				             &zu[nx + 1][0], nx2 * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 72, MPI_COMM_WORLD, &status);
//
//				MPI_Sendrecv(&xwd[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 63,
//				             &xd[nx + 1][0], nx2 * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 73, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&ywd[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 64,
//				             &yd[nx + 1][0], nx2 * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 74, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&zwd[1][0], nx * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 65,
//				             &zd[nx + 1][0], nx2 * (nz + 1), MPI_DOUBLE, myid + xln[xlln], 75, MPI_COMM_WORLD, &status);
//			}
		} else
			nx2 = 0;
		/*z����ǰ���������������*/

		if(rm == 1 && yl == 0)
			jtt = 24;
		else
			jtt = ny;

		//�������븺���±������(xf, yf, zf, xb, yb, zb),�˴���temp1����Ϊ��̬����,ԭ��ΧΪ (-nxyp:2*nxyp,-jtt:jt0)
		double **xf_temp1, **yf_temp1, **zf_temp1, **xb_temp1, **yb_temp1, **zb_temp1;
		double *xf_temp2[jtt + jt0 + 1], *yf_temp2[jtt + jt0 + 1], *zf_temp2[jtt + jt0 + 1], *xb_temp2[jtt + jt0 + 1], *yb_temp2[jtt + jt0 + 1], *zb_temp2[jtt + jt0 + 1];

		//��������ռ�Ķ�ά����
		xf_temp1 = malloc2D<double>(3 * nxyp + 1, jtt + jt0 + 1);//k=1
		yf_temp1 = malloc2D<double>(3 * nxyp + 1, jtt + jt0 + 1);
		zf_temp1 = malloc2D<double>(3 * nxyp + 1, jtt + jt0 + 1);
		xb_temp1 = malloc2D<double>(3 * nxyp + 1, jtt + jt0 + 1);//k=nz
		yb_temp1 = malloc2D<double>(3 * nxyp + 1, jtt + jt0 + 1);
		zb_temp1 = malloc2D<double>(3 * nxyp + 1, jtt + jt0 + 1);
		//�����ά��
		for(i = 0; i < jtt + jt0 + 1; i++) {
			xf_temp2[i] = xf_temp1[i] + nxyp;
			yf_temp2[i] = yf_temp1[i] + nxyp;
			zf_temp2[i] = zf_temp1[i] + nxyp;
			xb_temp2[i] = xb_temp1[i] + nxyp;
			yb_temp2[i] = yb_temp1[i] + nxyp;
			zb_temp2[i] = zb_temp1[i] + nxyp;
		}
		//�����ά��
		xf = xf_temp2 + jtt;
		yf = yf_temp2 + jtt;
		zf = zf_temp2 + jtt;
		xb = xb_temp2 + jtt;
		yb = yb_temp2 + jtt;
		zb = zb_temp2 + jtt;

		if(xlln == 1) {
			//�������븺���±������(xwf, ywf, zwf, xwb, ywb, zwb),�˴���temp1����Ϊ��̬���� ,ԭ��ΧΪ(1:nxyp,-jtt:jt0)
			double **xwf_temp1, **ywf_temp1, **zwf_temp1, **xwb_temp1, **ywb_temp1, **zwb_temp1;
			double *xwf_temp2[jtt + jt0 + 1], *ywf_temp2[jtt + jt0 + 1], *zwf_temp2[jtt + jt0 + 1], *xwb_temp2[jtt + jt0 + 1], *ywb_temp2[jtt + jt0 + 1], *zwb_temp2[jtt + jt0 + 1];

			//��������ռ�Ķ�ά����
			xwf_temp1 = malloc2D<double>(nxyp + 1, jtt + jt0 + 1);//k=1
			ywf_temp1 = malloc2D<double>(nxyp + 1, jtt + jt0 + 1);
			zwf_temp1 = malloc2D<double>(nxyp + 1, jtt + jt0 + 1);
			xwb_temp1 = malloc2D<double>(nxyp + 1, jtt + jt0 + 1);//k=nz
			ywb_temp1 = malloc2D<double>(nxyp + 1, jtt + jt0 + 1);
			zwb_temp1 = malloc2D<double>(nxyp + 1, jtt + jt0 + 1);
			//�����ά��
			for(i = 0; i < jtt + jt0 + 1; i++) {
				xwf_temp2[i] = xwf_temp1[i];
				ywf_temp2[i] = ywf_temp1[i];
				zwf_temp2[i] = zwf_temp1[i];
				xwb_temp2[i] = xwb_temp1[i];
				ywb_temp2[i] = ywb_temp1[i];
				zwb_temp2[i] = zwb_temp1[i];
			}
			//�����ά��
			xwf = xwf_temp2 + jtt;
			ywf = ywf_temp2 + jtt;
			zwf = zwf_temp2 + jtt;
			xwb = xwb_temp2 + jtt;
			ywb = ywb_temp2 + jtt;
			zwb = zwb_temp2 + jtt;

			for(j = jb; j < jt + 1; j++) { //Ĭ��z���򲻷ֿ�
				for(i = 1; i < nxyp + 1; i++) {
					xwf[j][i] = 0.25 * (x[1][j][i] + x[1][j][i + 1] + x[1][j + 1][i] + x[1][j + 1][i + 1]);    //ǰfront��������1������
					ywf[j][i] = 0.25 * (y[1][j][i] + y[1][j][i + 1] + y[1][j + 1][i] + y[1][j + 1][i + 1]);
					zwf[j][i] = 0.25 * (z[1][j][i] + z[1][j][i + 1] + z[1][j + 1][i] + z[1][j + 1][i + 1]);
					xwb[j][i] = 0.25 * (x[nz + 1][j][i] + x[nz + 1][j][i + 1] + x[nz + 1][j + 1][i] + x[nz + 1][j + 1][i + 1]);    //��back��������nz+1������
					ywb[j][i] = 0.25 * (y[nz + 1][j][i] + y[nz + 1][j][i + 1] + y[nz + 1][j + 1][i] + y[nz + 1][j + 1][i + 1]);
					zwb[j][i] = 0.25 * (z[nz + 1][j][i] + z[nz + 1][j][i + 1] + z[nz + 1][j + 1][i] + z[nz + 1][j + 1][i + 1]);
				}
			}
//		/*�̱�����ҶƬ���´��ݹ̱�ֵ*/
//		//�˴���Ϊ����������ȡѭ����ʽ����
			if(yl == 0) {
//			for(int iter = 1; iter < nxyp + 1; iter++) {
//				MPI_Sendrecv(&xwf[iter][jb], jt - jb + 1, MPI_DOUBLE, myid + xln[xlln], 360,
//				             &xwf[iter][jt + jb], jtt - jb + 1, MPI_DOUBLE, myid + xln[xlln], 460, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&ywf[iter][jb], jt - jb + 1, MPI_DOUBLE, myid + xln[xlln], 370,
//				             &ywf[iter][jt + jb], jtt - jb + 1, MPI_DOUBLE, myid + xln[xlln], 470, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&zwf[iter][jb], jt - jb + 1, MPI_DOUBLE, myid + xln[xlln], 380,
//				             &zwf[iter][jt + jb], jtt - jb + 1, MPI_DOUBLE, myid + xln[xlln], 480, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&xwb[iter][jb], jt - jb + 1, MPI_DOUBLE, myid + xln[xlln], 361,
//				             &xwb[iter][jt + jb], jtt - jb + 1, MPI_DOUBLE, myid + xln[xlln], 461, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&ywb[iter][jb], jt - jb + 1, MPI_DOUBLE, myid + xln[xlln], 371,
//				             &ywb[iter][jt + jb], jtt - jb + 1, MPI_DOUBLE, myid + xln[xlln], 471, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&zwb[iter][jb], jt - jb + 1, MPI_DOUBLE, myid + xln[xlln], 381,
//				             &zwb[iter][jt + jb], jtt - jb + 1, MPI_DOUBLE, myid + xln[xlln], 481, MPI_COMM_WORLD, &status);
//			}
//
				ny3 = jb;
				ny4 = jt0;
			} else if(yl == 1) {
//			for(int iter = 1; iter < nxyp + 1; iter++) {
//				MPI_Sendrecv(&xwf[iter][jb], jt - jb + 1, MPI_DOUBLE, myid - xln[xlln], 460,
//				             &xwf[iter][-jtt], jtt - jb + 1, MPI_DOUBLE, myid - xln[xlln], 360, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&ywf[iter][jb], jt - jb + 1, MPI_DOUBLE, myid - xln[xlln], 470,
//				             &ywf[iter][-jtt], jtt - jb + 1, MPI_DOUBLE, myid - xln[xlln], 370, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&zwf[iter][jb], jt - jb + 1, MPI_DOUBLE, myid - xln[xlln], 480,
//				             &zwf[iter][-jtt], jtt - jb + 1, MPI_DOUBLE, myid - xln[xlln], 380, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&xwb[iter][jb], jt - jb + 1, MPI_DOUBLE, myid - xln[xlln], 461,
//				             &xwb[iter][-jtt], jtt - jb + 1, MPI_DOUBLE, myid - xln[xlln], 361, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&ywb[iter][jb], jt - jb + 1, MPI_DOUBLE, myid - xln[xlln], 471,
//				             &ywb[iter][-jtt], jtt - jb + 1, MPI_DOUBLE, myid - xln[xlln], 371, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&zwb[iter][jb], jt - jb + 1, MPI_DOUBLE, myid - xln[xlln], 481,
//				             &zwb[iter][-jtt], jtt - jb + 1, MPI_DOUBLE, myid - xln[xlln], 381, MPI_COMM_WORLD, &status);
//			}
				ny3 = - jtt;
				ny4 = jt;
			}

			for(j = ny3; j < ny4 + 1; j++) {
				for(i = 1; i < nxyp + 1; i++) {
					xb[j][i] = xwb[j][i];
					yb[j][i] = ywb[j][i];
					zb[j][i] = zwb[j][i];
					xf[j][i] = xwf[j][i];
					yf[j][i] = ywf[j][i];
					zf[j][i] = zwf[j][i];
				}
			}
//
//		/*�����������͹̱������ݹ̱�ֵ*/
			if(xl == 0) {
				if(yl == 0)
					j = 0;
				else
					j = 1;
//				for(i = 1; i < xln[0] + 1; i++) {
//					MPI_SEND(ny3, 1, MPI_INT, myid - xln[j] - i, 210, MPI_COMM_WORLD);
//					MPI_SEND(ny4, 1, MPI_INT, myid - xln[j] - i, 220, MPI_COMM_WORLD);
//					//�˴���Ϊ����������ȡѭ����ʽ����
//					for(int iter = 1; iter < nxyp + 1; iter++) {
//						MPI_SEND(&xwb[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - xln[j] - i, 230, MPI_COMM_WORLD);
//						MPI_SEND(&ywb[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - xln[j] - i, 240, MPI_COMM_WORLD);
//						MPI_SEND(&zwb[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - xln[j] - i, 250, MPI_COMM_WORLD);
//						MPI_SEND(&xwf[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - xln[j] - i, 260, MPI_COMM_WORLD);
//						MPI_SEND(&ywf[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - xln[j] - i, 270, MPI_COMM_WORLD);
//						MPI_SEND(&zwf[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - xln[j] - i, 280, MPI_COMM_WORLD);
//					}
//
//				}
//				for(int iter = 1; iter < nxyp + 1; iter++) {
//					MPI_Sendrecv(&xwb[iter][ny3], nxyp * (ny4 - ny3 + 1), MPI_DOUBLE, myid + 1, 231,
//					             &xb[nxyp + iter][ny3], nxyp * (ny4 - ny3 + 1), MPI_DOUBLE, myid + 1, 241, MPI_COMM_WORLD, &status);
//					MPI_Sendrecv(&ywb[iter][ny3], nxyp * (ny4 - ny3 + 1), MPI_DOUBLE, myid + 1, 232,
//					             &yb[nxyp + iter][ny3], nxyp * (ny4 - ny3 + 1), MPI_DOUBLE, myid + 1, 242, MPI_COMM_WORLD, &status);
//					MPI_Sendrecv(&zwb[iter][ny3], nxyp * (ny4 - ny3 + 1), MPI_DOUBLE, myid + 1, 233,
//					             &zb[nxyp + iter][ny3], nxyp * (ny4 - ny3 + 1), MPI_DOUBLE, myid + 1, 243, MPI_COMM_WORLD, &status);
//					MPI_Sendrecv(&xwf[iter][ny3], nxyp * (ny4 - ny3 + 1), MPI_DOUBLE, myid + 1, 234,
//					             &xf[nxyp + iter][ny3], nxyp * (ny4 - ny3 + 1), MPI_DOUBLE, myid + 1, 244, MPI_COMM_WORLD, &status);
//					MPI_Sendrecv(&ywf[iter][ny3], nxyp * (ny4 - ny3 + 1), MPI_DOUBLE, myid + 1, 235,
//					             &yf[nxyp + iter][ny3], nxyp * (ny4 - ny3 + 1), MPI_DOUBLE, myid + 1, 245, MPI_COMM_WORLD, &status);
//					MPI_Sendrecv(&zwf[iter][ny3], nxyp * (ny4 - ny3 + 1), MPI_DOUBLE, myid + 1, 236,
//					             &zf[nxyp + iter][ny3], nxyp * (ny4 - ny3 + 1), MPI_DOUBLE, myid + 1, 246, MPI_COMM_WORLD, &status);
//				}
//
				nx3 = 1;
				nx4 = 2 * nxyp;
			}
			if(xl == xln[1] - 1) {
				if(yl == 0)
					j = 1;
				else
					j = 2;
//				for(i = 1; i < xln[2] + 1; i++) {
//					MPI_SEND(ny3, 1, MPI_INT, myid + xln[j] + i, 210, MPI_COMM_WORLD);
//					MPI_SEND(ny4, 1, MPI_INT, myid + xln[j] + i, 220, MPI_COMM_WORLD);
//					for(int iter = 1; iter < nxyp + 1; iter++) {
//						MPI_SEND(&xwb[iter][ny3], ny4 - ny3 + 1), MPI_DOUBLE, myid + xln[j] + i, 230, MPI_COMM_WORLD);
//						MPI_SEND(&ywb[iter][ny3], ny4 - ny3 + 1), MPI_DOUBLE, myid + xln[j] + i, 240, MPI_COMM_WORLD);
//						MPI_SEND(&zwb[iter][ny3], ny4 - ny3 + 1), MPI_DOUBLE, myid + xln[j] + i, 250, MPI_COMM_WORLD);
//						MPI_SEND(&xwf[iter][ny3], ny4 - ny3 + 1), MPI_DOUBLE, myid + xln[j] + i, 260, MPI_COMM_WORLD);
//						MPI_SEND(&ywf[iter][ny3], ny4 - ny3 + 1), MPI_DOUBLE, myid + xln[j] + i, 270, MPI_COMM_WORLD);
//						MPI_SEND(&zwf[iter][ny3], ny4 - ny3 + 1), MPI_DOUBLE, myid + xln[j] + i, 280, MPI_COMM_WORLD);
//					}
//				}
//				for(int iter = 0; iter < nxyp + 1; iter++) {
//					MPI_Sendrecv(&xwb[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - 1, 241,
//					             &xb[-nxyp + iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - 1, 231, MPI_COMM_WORLD, &status);
//					MPI_Sendrecv(&ywb[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - 1, 242,
//					             &yb[-nxyp + iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - 1, 232, MPI_COMM_WORLD, &status);
//					MPI_Sendrecv(&zwb[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - 1, 243,
//					             &zb[-nxyp + iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - 1, 233, MPI_COMM_WORLD, &status);
//					MPI_Sendrecv(&xwf[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - 1, 244,
//					             &xf[-nxyp + iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - 1, 234, MPI_COMM_WORLD, &status);
//					MPI_Sendrecv(&ywf[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - 1, 245,
//					             &yf[-nxyp + iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - 1, 235, MPI_COMM_WORLD, &status);
//					MPI_Sendrecv(&zwf[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - 1, 246,
//					             &zf[-nxyp + iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, myid - 1, 236, MPI_COMM_WORLD, &status);
//				}
//
				nx3 = - nxyp;
				nx4 = nxyp;
			} else {
				if(xlln == 0) {
					if(yl == 0)
						j = xln[0] * yln + lbm * numpp;
					else if(yl == 1)
						j = xln[0] * yln + xln[1] + lbm * numpp;
				} else if(xlln == 2) {
					if(yl == 0)
						j = xln[0] * yln + xln[1] - 1 + lbm * numpp;
					else if(yl == 1)
						j = xln[0] * yln + xln[1] * yln - 1 + lbm * numpp;
				}
			}
//			MPI_RECV(ny3, 1, MPI_INT, j, 210, MPI_COMM_WORLD, &status);
//			MPI_RECV(ny4, 1, MPI_INT, j, 220, MPI_COMM_WORLD, &status);
//			for(int iter = 0 ; iter < nxyp + 1; iter++) {
//				MPI_RECV(&xb[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, j, 230, MPI_COMM_WORLD, &status);
//				MPI_RECV(&yb[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, j, 240, MPI_COMM_WORLD, &status);
//				MPI_RECV(&zb[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, j, 250, MPI_COMM_WORLD, &status);
//				MPI_RECV(&xf[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, j, 260, MPI_COMM_WORLD, &status);
//				MPI_RECV(&yf[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, j, 270, MPI_COMM_WORLD, &status);
//				MPI_RECV(&zf[iter][ny3], ny4 - ny3 + 1, MPI_DOUBLE, j, 280, MPI_COMM_WORLD, &status);
//			}
//
//		}
//
//		/*�󵽱�����̾���*/
			for (int k = 1; k<nz + 1; k++)
				for (int j = 1; j<ny + 1; j++)
					for (int i = 1; i<nx + 1; i++) {
						dy = pow(10.0, 20);
						//��������¹̱ھ���,��ȷ��ĳ�ռ���¶�Ӧ�Ḻ́ڵ����귶Χ
						iil = i - 10;
						iil = max(iil, -nx1);
						iir = iil + 20;
						iir = min(iir, nx + nx2);

						kkl = k - 10;
						kkl = max(kkl, 1);
						kkr = kkl + 20;
						kkr = min(kkr, nz);
						for (kk = kkl; kk < kkr + 1; kk++)
							for (ii = iil; ii < iir + 1; ii++) {
								if(ii == 0) continue;

								dy1 = pow(xx00[k][j][i] - xd[kk][ii], 2) + pow(yy00[k][j][i] - yd[kk][ii], 2) + pow(zz00[k][j][i] - zd[kk][ii], 2);
								dy2 = pow(xx00[k][j][i] - xu[kk][ii], 2) + pow(yy00[k][j][i] - yu[kk][ii], 2) + pow(zz00[k][j][i] - zu[kk][ii], 2);

								d = min(dy1, dy2);
								if (d < dy)
									dy = d;
							}
						/*�����z�����ǰ��̱ھ���,ͬ�Ϸ�����ֻ��Ҫ������*/
						dz = pow(10.0, 20);

						if(xlln == 1) {
							iib = i - 10;
							iib = max(iib, nx3);
							iit = iib + 20;
							iit = min(iit, nx4);
						} else if(xlln == 0) {
							iib = 1;
							iit = iib + 20;
						} else if(xlln == 2) {
							iib = nxyp - 20;
							iit = nxyp;
						}

						if (j < jb) {
							jjb = jb;
							jjt = jb + 20;
						} else if ((j >= jb) && (j <= jt)) {
							jjb = j - 10;
							jjb = max(jjb, ny3);
							jjt = jjb + 20;
							jjt = min(jjt, ny4);
						} else {
							jjb = jt - 20;
							jjt = jt;
						}

						for (jj = jjb; jj < jjt + 1; jj++)
							for (ii = iib; ii<iit + 1; ii++) {
								if(ii == 0 || jj == 0) continue;

								dz1 = (pow((xx00[k][j][i] - xf[jj][ii]), 2) + pow((yy00[k][j][i] - yf[jj][ii]), 2) + pow((zz00[k][j][i] - zf[jj][ii]), 2));
								dz2 = (pow((xx00[k][j][i] - xb[jj][ii]), 2) + pow((yy00[k][j][i] - yb[jj][ii]), 2) + pow((zz00[k][j][i] - zb[jj][ii]), 2));

								d = min(dz1, dz2);
								if (d < dz)
									dz = d;
							}
						dmini[k][j][i] = sqrt(min(dy, dz));
					}
		}
	}

	/*
	 */
	void geo() {
		double temp, x24, y24, z24, x31, y31, z31;

		for (k = 1; k < nz + 2; k++)
			for (j = 1; j<ny + 1; j++)
				for (i = 1; i<nx + 1; i++) {
					x24 = x[k][j + 1][i] - x[k][j][i + 1];
					y24 = y[k][j + 1][i] - y[k][j][i + 1];
					z24 = z[k][j + 1][i] - z[k][j][i + 1];

					x31 = x[k][j + 1][i + 1] - x[k][j][i];
					y31 = y[k][j + 1][i + 1] - y[k][j][i];
					z31 = z[k][j + 1][i + 1] - z[k][j][i];

					s1x[k][j][i] = 0.5 * (y24 * z31 - z24 * y31);
					s1y[k][j][i] = 0.5 * (z24 * x31 - x24 * z31);
					s1z[k][j][i] = 0.5 * (x24 * y31 - y24 * x31);

					xx01[k][j][i] = 0.25 * (x[k][j][i] + x[k][j + 1][i] + x[k][j][i + 1] + x[k][j + 1][i + 1]);
					yy01[k][j][i] = 0.25 * (y[k][j][i] + y[k][j + 1][i] + y[k][j][i + 1] + y[k][j + 1][i + 1]);
					zz01[k][j][i] = 0.25 * (z[k][j][i] + z[k][j + 1][i] + z[k][j][i + 1] + z[k][j + 1][i + 1]);
				}

		for (k = 1; k<nz + 1; k++)
			for (j = 1; j<ny + 1; j++)
				for (i = 1; i<nx + 2; i++) {
					x24 = x[k + 1][j + 1][i] - x[k][j][i];
					y24 = y[k + 1][j + 1][i] - y[k][j][i];
					z24 = z[k + 1][j + 1][i] - z[k][j][i];

					x31 = x[k][j + 1][i] - x[k + 1][j][i];
					y31 = y[k][j + 1][i] - y[k + 1][j][i];
					z31 = z[k][j + 1][i] - z[k + 1][j][i];

					s2x[k][j][i] = 0.5 * (y24 * z31 - z24 * y31);
					s2y[k][j][i] = 0.5 * (z24 * x31 - x24 * z31);
					s2z[k][j][i] = 0.5 * (x24 * y31 - y24 * x31);

					xx02[k][j][i] = 0.25 * (x[k][j][i] + x[k][j + 1][i] + x[k + 1][j][i] + x[k + 1][j + 1][i]);
					yy02[k][j][i] = 0.25 * (y[k][j][i] + y[k][j + 1][i] + y[k + 1][j][i] + y[k + 1][j + 1][i]);
					zz02[k][j][i] = 0.25 * (z[k][j][i] + z[k][j + 1][i] + z[k + 1][j][i] + z[k + 1][j + 1][i]);
				}

		for (k = 1; k < nz + 1; k++)
			for (j = 1; j<ny + 2; j++)
				for (i = 1; i<nx + 1; i++) {

					x24 = x0[k + 1][j][i + 1] - x[k][j][i];
					x31 = x[k + 1][j][i] - x[k][j][i + 1];

					y24 = y[k + 1][j][i + 1] - y[k][j][i];
					y31 = y[k + 1][j][i] - y[k][j][i + 1];

					z24 = z[k + 1][j][i + 1] - z[k][j][i];
					z31 = z[k + 1][j][i] - z[k][j][i + 1];

					s3x[k][j][i] = 0.5 * (y24 * z31 - z24 * y31);
					s3y[k][j][i] = 0.5 * (z24 * x31 - x24 * z31);
					s3z[k][j][i] = 0.5 * (x24 * y31 - y24 * x31);

					xx03[k][j][i] = 0.25 * (x[k][j][i] + x[k][j][i + 1] + x[k + 1][j][i] + x[k + 1][j][i + 1]);
					yy03[k][j][i] = 0.25 * (y[k][j][i] + y[k][j][i + 1] + y[k + 1][j][i] + y[k + 1][j][i + 1]);
					zz03[k][j][i] = 0.25 * (z[k][j][i] + z[k][j][i + 1] + z[k + 1][j][i] + z[k + 1][j][i + 1]);

				}

		for (k = 1; k < nz + 1; k++)
			for (j = 1; j<ny + 1; j++)
				for (i = 1; i<nx + 1; i++) {

					x24 = x[k + 1][j + 1][i + 1] - x[k][j][i];
					y24 = y[k + 1][j + 1][i + 1] - y[k][j][i];
					z24 = z[k + 1][j + 1][i + 1] - z[k][j][i];

					vv0[k][j][i] = -(x24 * (s1x[k][j][i] + s2x[k][j][i] + s3x[k][j][i]) +
					                 y24 * (s1y[k][j][i] + s2y[k][j][i] + s3y[k][j][i]) +
					                 z24 * (s1z[k][j][i] + s2z[k][j][i] + s3z[k][j][i])) / 3.0;

					if(vv0[k][j][i] < 0) {
						cout<<"vv0 <= 0" << i << j << k << endl;
					}

					xx0[k][j][i] = 0.125 * (x[k][j][i] + x[k][j][i + 1] + x[k + 1][j][i] + x[k + 1][j][i + 1] +
					                        x[k][j + 1][i] + x[k][j + 1][i + 1] + x[k + 1][j + 1][i] + x[k + 1][j + 1][i + 1]);

					yy0[k][j][i] = 0.125 * (y[k][j][i] + y[k][j][i + 1] + y[k + 1][j][i] + y[k + 1][j][i + 1] +
					                        y[k][j + 1][i] + y[k][j + 1][i + 1] + y[k + 1][j + 1][i] + y[k + 1][j + 1][i + 1]);

					zz0[k][j][i] = 0.125 * (z[k][j][i] + z[k][j][i + 1] + z[k + 1][j][i] + z[k + 1][j][i + 1] +
					                        z[k][j + 1][i] + z[k][j + 1][i + 1] + z[k + 1][j + 1][i] + z[k + 1][j + 1][i + 1]);
				}

		for (k = 1; k<nz + 1; k++)
			for (j = 1; j<ny + 1; j++) {
				s2x[k][j][0] = s2x[k][j][1];
				s2y[k][j][0] = s2y[k][j][1];
				s2z[k][j][0] = s2z[k][j][1];

				s2x[k][j][nx + 2] = s2x[k][j][nx + 1];
				s2y[k][j][nx + 2] = s2y[k][j][nx + 1];
				s2z[k][j][nx + 2] = s2z[k][j][nx + 1];
			}

		for (k = 1; k<nz + 1; k++)
			for (i = 1; i<nx + 1; i++) {
				s3x[k][0][i] = s3x[k][1][i];
				s3y[k][0][i] = s3y[k][1][i];
				s3z[k][0][i] = s3z[k][1][i];

				s3x[k][ny + 2][i] = s3x[k][ny + 1][i];
				s3y[k][ny + 2][i] = s3y[k][ny + 1][i];
				s3z[k][ny + 2][i] = s3z[k][ny + 1][i];
			}

		for (j = 1; j<ny + 1; j++)
			for (i = 1; i<nx + 1; i++) {
				s1x[0][j][i] = s1x[1][j][i];
				s1y[0][j][i] = s1y[1][j][i];
				s1z[0][j][i] = s1z[1][j][i];

				s1x[nz + 2][j][i] = s1x[nz + 1][j][i];
				s1y[nz + 2][j][i] = s1y[nz + 1][j][i];
				s1z[nz + 2][j][i] = s1z[nz + 1][j][i];
			}

		/*�������񶥵�����ֵ,���������һ���������*/
		if(slidm == 1) {
			for(k = 1; k < nz + 2; k++) {
				for(j = 1; j < ny + 2; j++) {
					x[k][j][nx + 2] = 2 * x[k][j][nx + 1] - x[k][j][nx];//���񶥵�����
					y[k][j][nx + 2] = 2 * y[k][j][nx + 1] - y[k][j][nx];
					z[k][j][nx + 2] = 2 * z[k][j][nx + 1] - z[k][j][nx];
				}
			}
			for(k = 1; k < nz + 1; k++)
				for(j = 1; j < ny + 1; j++)
					vv[k][j] = vv0[k][j][nx];

		} else if(slidm == 2) {
			for(k = 1; k < nz + 2; k++) {
				for(j = 1; j < ny + 2; j++) {
					x[k][j][0] = 2 * x[k][j][1] - x[k][j][2];//���񶥵�����
					y[k][j][0] = 2 * y[k][j][1] - y[k][j][2];
					z[k][j][0] = 2 * z[k][j][1] - z[k][j][2];
				}
			}
			for(k = 1; k < nz + 1; k++)
				for(j = 1; j < ny + 1; j++)
					vv[k][j] = vv0[k][j][1];
		}
	}

///* �������� 1
// * ����y���������Ϊ������������,�����ǰ������Ҷ�߰ٷֱȶ�Ӧ��ʵ�ʾ���ֵ
// */
	void spaa() {
		double t1, t2, temp;
		double *sss, *hx, *hy, *hz, *hr;
		double ***xff, ***yff, ***zff;

		ny0 = ny * yln;

		xff = malloc3D<double>(nx + 2, ny0 + 2, nz + 2);//�����yln=1�µ�����
		yff = malloc3D<double>(nx + 2, ny0 + 2, nz + 2);
		zff = malloc3D<double>(nx + 2, ny0 + 2, nz + 2);

		if(rm == 1) {
			nx1 = nxyp;//��ǰ����ǰ��ÿ����������
			for(k = 1; k < nz + 2; k++) {
				for(j = 1; j < ny0 + 2; j++) {
					for(i = 1; i < nx + 2; i++) {
						xff[k][j][i] = xf0[zl * nz + k][j][xll * nx1 + i];
						yff[k][j][i] = yf0[zl * nz + k][j][xll * nx1 + i];
						zff[k][j][i] = zf0[zl * nz + k][j][xll * nx1 + i];
					}
				}
			}
		} else {
			nx1 = 24;//��ǰ����ǰ��ÿ����������
			for(k = 1; k < nz + 2; k++) {
				for(j = 1; j < ny0 + 2; j++) {
					for(i = 1; i < nx + 2; i++) {
						xff[k][j][i] = xf0[zl * nz + k][j][xll * nx - nx + nx1 + i];
						yff[k][j][i] = yf0[zl * nz + k][j][xll * nx - nx + nx1 + i];
						zff[k][j][i] = zf0[zl * nz + k][j][xll * nx - nx + nx1 + i];
					}
				}
			}
		}

		temp = 2 / double(lb[rm]) * pi * lm;//��ת������
		cor = cos(temp);						//0             ^
		sir = sin(temp);						//1             ^        ҶƬ������ת
		for(k = 1; k < nz + 2; k++) 			//nz            ^
			for(j = 1; j < ny0 + 2; j++) 		//nz+1          1
				for(i = 1; i < nx + 2; i++) {
					t1 = yff[k][j][i];
					t2 = zff[k][j][i];
					yff[k][j][i] = t1 * cor + t2 * sir;//��ҶƬ��ת������ת,���෴����Ϊ����0��nz+1���꣬1-nz��˳����ҶƬ��ת�����෴
					zff[k][j][i] = - t1 * sir + t2 * cor;
				}



		Dmalloc1(sss,ny0 + 2);
		Dmalloc1(hx, ny0 + 2);
		Dmalloc1(hy, ny0 + 2);
		Dmalloc1(hz, ny0 + 2);
		Dmalloc1(hr, ny0 + 2);
		if(spa1 < 50)//ֻ�Ǹ������жϱ�׼
			l = 0;
		else
			l = 1;
		if(spa2 < 50)//ֻ�Ǹ������жϱ�׼
			kk = 0;
		else
			kk = 1;

		for(k = 1; j < nz + 1; k++) {
			for(i = 1; i < nx + 1; i++) {
				sss[1] = 0;
				for(j = 1; j < ny0 + 1; j++) {
					hx[j] = xff[k][j][i];
					hy[j] = yff[k][j][i];
					hz[j] = zff[k][j][i];
					hr[j] = sqrt(hy[j] * hy[j] + hz[j] * hz[j]);
					if(j > 1)
						sss[j] = sss[j - 1] + sqrt( pow(hx[j] - hx[j - 1], 2) + pow(hr[j] - hr[j - 1], 2) );//����Ϊ���������������������
				}
				if(l==0 && yl==l) {
					t1 = sss[ny0 + 1] * spa1 / 100;
//				wl(sss, hx, t1, hxx1[k][i], ny0 + 1, 1);
//				wl(sss, hy, t1, hyy1[k][i], ny0 + 1, 1);
//				wl(sss, hz, t1, hzz1[k][i], ny0 + 1, 1);
				} else if (l==1 && yl==l) {
					t1 = sss[ny0 + 1] * spa1 / 100;
//				wl(sss, hx, t1, hxx1[k][i], ny0 + 1, 1);
//				wl(sss, hy, t1, hyy1[k][i], ny0 + 1, 1);
//				wl(sss, hz, t1, hzz1[k][i], ny0 + 1, 1);
				}
				if(kk==0 && yl==kk) {
					t2 = sss[ny0 + 1] * spa2 / 100;
//				wl(sss, hx, t2, hxx2[k][i], ny0 + 1, 1);
//				wl(sss, hy, t2, hyy2[k][i], ny0 + 1, 1);
//				wl(sss, hz, t2, hzz2[k][i], ny0 + 1, 1);
				} else if(kk==1 && yl==kk) {
					t2 = sss[ny0 + 1] * spa2 / 100;
//				wl(sss, hx, t2, hxx2[k][i], ny0 + 1, 1);
//				wl(sss, hy, t2, hyy2[k][i], ny0 + 1, 1);
//				wl(sss, hz, t2, hzz2[k][i], ny0 + 1, 1);
				}
			}
		}

		for(k = 1; k < nz + 1; k++) {
			for(i = 1; i < nx + 1; i++) {
				sss[1] = 0;
				for(j = 1; j < ny0 + 1; j++) {
					hx[j] = 0.125 * (xff[k][j][i] + xff[k][j][i + 1] + xff[k + 1][j][i] + xff[k + 1][j][i + 1] +
					                 xff[k][j + 1][i] + xff[k][j + 1][i + 1] + xff[k + 1][j + 1][i] + xff[k + 1][j + 1][i + 1]);
					y1 = 0.125 * (yff[k][j][i] + yff[k][j][i + 1] + yff[k + 1][j][i] + yff[k + 1][j][i + 1] +
					              yff[k][j + 1][i] + yff[k][j + 1][i + 1] + yff[k + 1][j + 1][i] + yff[k + 1][j + 1][i + 1]);
					z1 = 0.125 * (zff[k][j][i] + zff[k][j][i + 1] + zff[k + 1][j][i] + zff[k + 1][j][i + 1] +
					              zff[k][j + 1][i] + zff[k][j + 1][i + 1] + zff[k + 1][j + 1][i] + zff[k + 1][j + 1][i + 1]);
					hr[j] = sqrt(y1 * y1 + z1 * z1);
					if(j > 1)
						sss[j] = sss[j - 1] + sqrt( pow(hx[j] - hx[j - 1], 2) + pow(hr[j] - hr[j - 1], 2) );//����Ϊ���������������������
				}
				hspa[k][i] = sss[ny + 1];
				if(l == 0 && yl == l)
					spa01[k][i] = sss[ny0] * spa1 / 100;
				else if (l == 1 && yl == l)
					spa01[k][i] = sss[ny0] * spa1 / 100;
				if(kk==0 && yl == kk)
					spa02[k][i] = sss[ny0] * spa2 / 100;
				else if (kk == 1 && yl == kk)
					spa02[k][i] = sss[ny0] * spa2 / 100;
			}
		}
		/*Ϊ����������Ҷ�߷ֲ������*/
		if(rm == 2 && xll == 4 && lm == 0) {
			Dmalloc1(temp0, ny0 + 1);
			i = nx;
			k = 1;
			for(j = 1; j < ny0 + 1; j++) {
				y1 = 0.125 * (yff[k][j][i] + yff[k][j][i + 1] + yff[k + 1][j][i] + yff[k + 1][j][i + 1] +
				              yff[k][j + 1][i] + yff[k][j + 1][i + 1] + yff[k + 1][j + 1][i] + yff[k + 1][j + 1][i + 1]);
				z1 = 0.125 * (zff[k][j][i] + zff[k][j][i + 1] + zff[k + 1][j][i] + zff[k + 1][j][i + 1] +
				              zff[k][j + 1][i] + zff[k][j + 1][i + 1] + zff[k + 1][j + 1][i] + zff[k + 1][j + 1][i + 1]);
				hr[j] = sqrt(y1 * y1 + z1 * z1);
			}
			for(j = 1; j < ny0 + 1; j++) {
				temp0[j] = (hr[j] - hr[1]) / (hr[ny0] - hr[1]) * 100;
			}
		}
	}

	void slid() {
		double temp;
		double ***ys, ***zs;

		temp = delt * double(n - 1 + (m - 1) * nt) * rpm[1];//��ǰ�⡢��ѭ������ת����ת���ĽǶȡ�rpmΪ��ֵ���ʺ��涯ҶҪȡ����
		if(rm == 1) {
			i = nx - 1;
			temp = - temp;//���ζ�Ҷ�����񶥵�ӳ�䵽��Ҷ�ϣ�����Ҫ�������->��������a,��ҶƬ��ת����ת
		} else if(rm == 2) {
			i = 1;
			temp = temp;  //���ξ�Ҷ�����񶥵�ӳ�䵽��Ҷ�ϣ�����Ҫ��������->���������r,��ҶƬ��ת����ת
		}

		ys = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		zs = malloc3D<double>(nx + 2, ny + 2, nz + 2);
		cor = cos(temp);//n��2��ʼ���ʸտ�ʼ����ʱ�Ѿ��ƶ���0.05��ͨ��
		sir = sin(temp);
		for(l = i; l < i + 3; l++)
			for(k = 1; k < nz + 2; k++)
				for(j = 1; j < ny + 2; j++) {
					ys[k][j][i] = y[k][j][i] * cor + z[k][j][i] * sir;
					zs[k][j][i] = - y[k][j][i] * sir + z[k][j][i] * cor;
				}

		//(4, ny + 2, nz + 2, lb[1] + lb[2])
		for(int iter = 0; iter < lb[1] + lb[2]; iter++)
		for(k = 1; k < nz + 2; k++)
			for(j = 1; j < ny + 2; j++)
				for(i = 1; i < 4; i++) {
					xxs[iter][k][j][i] = 0;
					yys[iter][k][j][i] = 0;
					zzs[iter][k][j][i] = 0;
			}
		//i�ķ���Ϊ���ά�ȣ���ά���ϲ�����������ֻ�ܲ���ѭ����ʽ���� 
		if(rm == 1) {
			for(j = yl + lb[rm] * numpp; j <= yl + lb[rm] * numpp + (lb[rm + 1] - 1) * numpp; j += numpp) {
				k = j / numpp;
				//ͨѶ��δ�޸ģ�����ȷ 
//				MPI_Sendrecv(&x[0][0][i], 4 * (ny + 2) * (nz + 2), MPI_DOUBLE, j, 5,
//				             &xxs[k][0][0][0], 4 * (ny + 2) * (nz + 2), MPI_DOUBLE, j, 15, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&ys[0][0][i], 4 * (ny + 2) * (nz + 2), MPI_DOUBLE, j, 6,
//				             &yys[k][0][0][0], 4 * (ny + 2) * (nz + 2), MPI_DOUBLE, j, 16, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&zs[0][0][i], 4 * (ny + 2) * (nz + 2), MPI_DOUBLE, j, 7,
//				             &zzs[k][0][0][0], 4 * (ny + 2) * (nz + 2), MPI_DOUBLE, j, 17, MPI_COMM_WORLD, &status);
			}
		} else if(rm == 2) {
			for(j = yl + 2 * xln[1] * yln; j <= yl + 2 * xln[1] * yln + (lb[rm - 1] - 1) * numpp; j += numpp) {
				k = j / numpp;
//				MPI_Sendrecv(&x[0][0][i], 4 * (ny + 2) * (nz + 2), MPI_DOUBLE, j, 15,
//				             &xxs[k][0][0][0], 4 * (ny + 2) * (nz + 2), MPI_DOUBLE, j, 5, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&ys[0][0][i], 4 * (ny + 2) * (nz + 2), MPI_DOUBLE, j, 16,
//				             &yys[k][0][0][0], 4 * (ny + 2) * (nz + 2), MPI_DOUBLE, j, 6, MPI_COMM_WORLD, &status);
//				MPI_Sendrecv(&zs[0][0][i], 4 * (ny + 2) * (nz + 2), MPI_DOUBLE, j, 17,
//				             &zzs[k][0][0][0], 4 * (ny + 2) * (nz + 2), MPI_DOUBLE, j, 7, MPI_COMM_WORLD, &status);
			}
		}
	}


}

//***********************************************
//***********************************************
//***********************************************
int main(int argc, char **argv) {
	cfd::init(argc, argv);

	cfd::ini();			//�����̶�����Ʋ�����ȫ����������
	cfd::allocation();	//�����ڴ�
	cfd::ini2();		//��1ʱ���ʼ��
	cfd::distance();	//��ϸ������ϣ�����������ĵ㵽�̱�����̵ľ���
	cfd::geo();			//�����β���
	cfd::spaa();		//����y���������Ϊ������������
	cfd::slid();		//�˺���ʵ���ϲ��ڴ˴������ڴ˴���Ϊ����֤�Ƿ���ȷ
	//cfd();			//spa1��spa2Ҷ�߽���������

	cfd::finished();
	return 0;
}
//***********************************************
//***********************************************
//***********************************************


