#define TESTMODE
#include "BezierCurveReconstruction.h"
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <fstream>
typedef Point2* BezierCurve;
#define MAXPOINTS	1000		/* The most points you can have */ // số điểm tối đa mà bạn có thể có


void DrawBezierCurve(int n, BezierCurve bcurve)
{
	int i;


	for (i = 0; i <= n; i++) {
		printf("%0.2f, %0.2f\n", bcurve[i].x, bcurve[i].y);
	}

}

/*
 *  B0, B1, B2, B3 :
 *	Bezier multipliers  //bội số Bezier
 */
static double B0(double	u)
{
	double tmp = 1.0 - u;
	return (tmp * tmp * tmp);
}


static double B1(double	u)
{
	double tmp = 1.0 - u;
	return (3 * u * (tmp * tmp));
}

static double B2(double	u)
{
	double tmp = 1.0 - u;
	return (3 * u * u * tmp);
}

static double B3(double	u)
{
	return (u * u * u);
}



static Vector2 V2AddII(Vector2 a, Vector2 b) {
	Vector2 c;
	c.x = a.x + b.x;  c.y = a.y + b.y;
	return (c);
}
static Vector2 V2ScaleIII(Vector2 v, double s) {
	Vector2 result;
	result.x = v.x * s; result.y = v.y * s;
	return (result);
}

static Vector2 V2SubII(Vector2 a, Vector2 b) {
	Vector2 c;
	c.x = a.x - b.x; c.y = a.y - b.y;
	return (c);
}

/* return the distance between two points */ // trả về khoảng cách giữa 2 điểm
double V2DistanceBetween2Points(Point2* a, Point2* b) {
	double dx = a->x - b->x;
	double dy = a->y - b->y;
	return(sqrt((dx * dx) + (dy * dy)));
}

/* negates the input vector and returns it */ // phủ định vector đầu vào và trả về nó
Vector2* V2Negate(Vector2* v) {
	v->x = -v->x;  v->y = -v->y;
	return(v);
}

double V2SquaredLength(Vector2* a) {
	return((a->x * a->x) + (a->y * a->y));
}


/* returns length of input vector */ //trả về độ dài của đầu vào vector
double V2Length(Vector2* a) {
	return(sqrt(V2SquaredLength(a)));
}

/* normalizes the input vector and returns it */ // bình thường hoá đầu vào vector và trả về nó
Vector2* V2Normalize(Vector2* v) {
	double len = V2Length(v);
	if (len != 0.0) { v->x /= len;  v->y /= len; }
	return(v);
}

/* scales the input vector to the new length and returns it */ // chia tỉ lệ vector đầu vào sang chiều dài mới và trả về nó
Vector2* V2Scale(Vector2* v, double newlen)
{
	double len = V2Length(v);
	if (len != 0.0)
	{
		v->x *= newlen / len;
		v->y *= newlen / len;
	}
	return (v);
}

/* return vector sum c = a+b */ // trả về vector tổng c = a+b
Vector2* V2Add(Vector2* a, Vector2* b, Vector2* c)
{
	c->x = a->x + b->x;
	c->y = a->y + b->y;
	return (c);
}

/* return vector difference c = a-b */ // trả về sự khác biệt vector c = a-b
Vector2* V2Sub(Vector2* a, Vector2* b, Vector2* c)
{
	c->x = a->x - b->x;
	c->y = a->y - b->y;
	return (c);
}

/* return the dot product of vectors a and b */ // trả về sản phẩm chấm của vectơ a và b
double V2Dot(Vector2* a, Vector2* b)
{
	return ((a->x * b->x) + (a->y * b->y));
}

/*
 *  ChordLengthParameterize : // độ dài hợp âm tham số hóa
 *  Assign parameter values to digitized points  // Gán các giá trị tham số cho các điểm được
 *  using relative distances between points. // số hóa bằng khoảng cách tương đối giữa các điểm
 */
static double* ChordLengthParameterize(Point2* d, int first, int last)
{
	int     i;
	double* u;         /*  Parameterization        */ // tham số hoá

	u = (double*)malloc((unsigned)(last - first + 1) * sizeof(double));

	u[0] = 0.0;
	for (i = first + 1; i <= last; i++) {
		u[i - first] = u[i - first - 1] +
			V2DistanceBetween2Points(&d[i], &d[i - 1]);
	}

	for (i = first + 1; i <= last; i++) {
		u[i - first] = u[i - first] / u[last - first];
	}

	return(u);
}



/*
 *  GenerateBezier : // tạo Bezier
 *  Use least-squares method to find Bezier control points for region.
 *  sử dụng phương pháp bình phương nhỏ nhất để tìm các điểm kiểm soát Bezier cho khu vực.
 */
static BezierCurve GenerateBezier(
	Point2* d, /*  Array of digitized points	*/                     // Mảng các điểm số hóa
	int first, int last, /*  Indices defining region	*/             // Chỉ số xác định vùng
	double* uPrime, /*  Parameter values for region */                 // Giá trị tham số cho vùng
	Vector2 tHat1, Vector2 tHat2) /*  Unit tangents at endpoints	*/ // Tiếp tuyến đơn vị tại các điểm cuối
{
	int i;
	Vector2 A[MAXPOINTS][2]; /* Precomputed rhs for eqn	*/ // ??
	int nPts; /* Number of pts in sub-curve */             // Số lượng điểm trong đường cong phụ
	double C[2][2]; /* Matrix C		*/                     // Ma trận C
	double X[2]; /* Matrix X			*/                 // Ma trận X
	double det_C0_C1, /* Determinants of matrices	*/     // Các yếu tố quyết định của ma trận
		det_C0_X,
		det_X_C1;
	double alpha_l, /* Alpha values, left and right	*/ // Giá trị Alpha, trái và phải
		alpha_r;
	Vector2 tmp; /* Utility variable		*/                 // Biến tiện ích
	BezierCurve bezCurve; /* RETURN bezier curve ctl pts	*/ // trả vở đường cong Bezier ctl pts

	bezCurve = (Point2*)malloc(4 * sizeof(Point2));
	nPts = last - first + 1;

	/* Compute the A's	*/ // Tính toán A's
	for (i = 0; i < nPts; i++)
	{
		Vector2 v1, v2;
		v1 = tHat1;
		v2 = tHat2;
		V2Scale(&v1, B1(uPrime[i]));
		V2Scale(&v2, B2(uPrime[i]));
		A[i][0] = v1;
		A[i][1] = v2;
	}

	/* Create the C and X matrices	*/ // Tạo ma trận C và X
	C[0][0] = 0.0;
	C[0][1] = 0.0;
	C[1][0] = 0.0;
	C[1][1] = 0.0;
	X[0] = 0.0;
	X[1] = 0.0;

	for (i = 0; i < nPts; i++)
	{
		C[0][0] += V2Dot(&A[i][0], &A[i][0]);
		C[0][1] += V2Dot(&A[i][0], &A[i][1]);
		/*					C[1][0] += V2Dot(&A[i][0], &A[i][1]);*/
		C[1][0] = C[0][1];
		C[1][1] += V2Dot(&A[i][1], &A[i][1]);

		tmp = V2SubII(d[first + i],
			V2AddII(
				V2ScaleIII(d[first], B0(uPrime[i])),
				V2AddII(
					V2ScaleIII(d[first], B1(uPrime[i])),
					V2AddII(
						V2ScaleIII(d[last], B2(uPrime[i])),
						V2ScaleIII(d[last], B3(uPrime[i]))))));

		X[0] += V2Dot(&A[i][0], &tmp);
		X[1] += V2Dot(&A[i][1], &tmp);
	}

	/* Compute the determinants of C and X	*/ // Tính các yếu tố quyết định của ma trận C và X
	det_C0_C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1];
	det_C0_X = C[0][0] * X[1] - C[0][1] * X[0];
	det_X_C1 = X[0] * C[1][1] - X[1] * C[0][1];

	/* Finally, derive alpha values	*/ // Cuối cùng, lấy giá trị alpha
	if (det_C0_C1 == 0.0)
	{
		det_C0_C1 = (C[0][0] * C[1][1]) * 10e-12;
	}
	alpha_l = det_X_C1 / det_C0_C1;
	alpha_r = det_C0_X / det_C0_C1;

	/*  If alpha negative, use the Wu/Barsky heuristic (see text) */
	// Nếu alpha âm, sử dụng Wu/Barsky heuristic
	/* (if alpha is 0, you get coincident control points that lead to
	* divide by zero in any subsequent NewtonRaphsonRootFind() call. */
	// Nếu alpha = 0, bạn nhận được các điểm kiểm soát trùng khớp dẫn đến
	// chia cho 0 trong bất kỳ lệnh gọi NewtonRaphsonRootFind () tiếp theo nào
	if (alpha_l < 1.0e-6 || alpha_r < 1.0e-6)
	{
		double dist = V2DistanceBetween2Points(&d[last], &d[first]) /
			3.0;

		bezCurve[0] = d[first];
		bezCurve[3] = d[last];
		V2Add(&bezCurve[0], V2Scale(&tHat1, dist), &bezCurve[1]);
		V2Add(&bezCurve[3], V2Scale(&tHat2, dist), &bezCurve[2]);
		return (bezCurve);
	}

	/*  First and last control points of the Bezier curve are */
	/*  positioned exactly at the first and last data points */
	/*  Control points 1 and 2 are positioned an alpha distance out */
	/*  on the tangent vectors, left and right, respectively */
	// Điểm kiểm soát đầu tiên và cuối cùng của đường cong Bezier được định vị
	// chính xác tại các điểm dữ liệu đầu tiên và cuối cùng
	//  Điểm kiểm soát 1 và 2 được định vị một khoảng cách alpha trên các vectơ tiếp tuyến, trái và phải, tương ứng
	bezCurve[0] = d[first];
	bezCurve[3] = d[last];
	V2Add(&bezCurve[0], V2Scale(&tHat1, alpha_l), &bezCurve[1]);
	V2Add(&bezCurve[3], V2Scale(&tHat2, alpha_r), &bezCurve[2]);
	return (bezCurve);
}

/*
 *  Bezier :
 *  	Evaluate a Bezier curve at a particular parameter value
 *      Đánh giá đường cong Bezier tại một giá trị tham số cụ thể
 */
static Point2 BezierII(
	int degree, /* The degree of the bezier curve	*/ //Độ của đường cong bezier
	Point2* V, /* Array of control points		*/     // Mảng các điểm kiểm soát
	double t) /* Parametric value to find point for	*/ //Giá trị tham số để tìm điểm cho
{
	int i, j;
	Point2 Q; /* Point on curve at parameter t	*/         //Điểm trên đường cong tại tham số t
	Point2* Vtemp; /* Local copy of control points		*/ //Bản sao các điểm kiểm soát

	/* Copy array	*/ // Sao chép mảng
	Vtemp = (Point2*)malloc((unsigned)((degree + 1) * sizeof(Point2)));
	for (i = 0; i <= degree; i++)
	{
		Vtemp[i] = V[i];
	}

	/* Triangle computation	*/ //Tam giác tính toán
	for (i = 1; i <= degree; i++)
	{
		for (j = 0; j <= degree - i; j++)
		{
			Vtemp[j].x = (1.0 - t) * Vtemp[j].x + t * Vtemp[j + 1].x;
			Vtemp[j].y = (1.0 - t) * Vtemp[j].y + t * Vtemp[j + 1].y;
		}
	}

	Q = Vtemp[0];
	free((void*)Vtemp);
	return Q;
}

/*
 * ComputeLeftTangent, ComputeRightTangent, ComputeCenterTangent : tính toán tiếp tuyến trái, tính toán tiếp tuyến phải
 tính toán tiếp tuyến giữa
 *Approximate unit tangents at endpoints and "center" of digitized curve
 Các tiếp tuyến đơn vị gần đúng tại các điểm cuối và "tâm" của đường cong số hóa
 */
static double ComputeMaxError(
	Point2* d,			/*  Array of digitized points	*/ // Mảng các điểm số hóa
	int		first, int last,		/*  Indices defining region	*/ // Chỉ số xác định vùng
	BezierCurve	bezCurve,		/*  Fitted Bezier curve		*/ //Đường cong Bezier được trang bị
	double* u,			/*  Parameterization of points	*/ //Tham số hóa điểm
	int* splitPoint)		/*  Point of maximum error	*/ // Điểm sai số tối đa
{
	int		i;
	double	maxDist;		/*  Maximum error		*/ // Lỗi tối đa
	double	dist;		/*  Current error		*/ // Lỗi hiện tại
	Point2	P;			/*  Point on curve		*/ // Điểm trên đường cong
	Vector2	v;			/*  Vector from point to curve	*/ // Vector từ điểm đến đường cong

	*splitPoint = (last - first + 1) / 2;
	maxDist = 0.0;
	for (i = first + 1; i < last; i++) {
		P = BezierII(3, bezCurve, u[i - first]);
		v = V2SubII(P, d[i]);
		dist = V2SquaredLength(&v);
		if (dist >= maxDist) {
			maxDist = dist;
			*splitPoint = i;
		}
	}
	return (maxDist);
}


/*
 *  NewtonRaphsonRootFind : công cụ newton-raph  cho việc tìm kiếm căn số
 *	Use Newton-Raphson iteration to find better root. : sử dụng công cụ đó lặp đi lặp lặp lại
	để tìm căn số tốt hơn
 */
static double NewtonRaphsonRootFind(
	BezierCurve	Q,			/*  Current fitted curve: đường cong thích hợp hiện tại	*/
	Point2 		P,		/*  Digitized point	: điểm số 2	*/
	double 		u)		/*  Parameter value for "P" : giá trị tham số cho P	*/
{
	double 		numerator, denominator;
	Point2 		Q1[3], Q2[2];	/*  Q' and Q'':	Q và Q'		*/
	Point2		Q_u, Q1_u, Q2_u; /*u evaluated at: u được xac định tại điểm Q, Q', & Q''	*/
	double 		uPrime;		/*  Improved u			*/
	int 		i;

	/* Compute Q(u)	: tính toán Q(u)*/
	Q_u = BezierII(3, Q, u);

	/* Generate control vertices for Q' : tạo các đỉnh điều khiển cho Q'	*/
	for (i = 0; i <= 2; i++) {
		Q1[i].x = (Q[i + 1].x - Q[i].x) * 3.0;
		Q1[i].y = (Q[i + 1].y - Q[i].y) * 3.0;
	}

	/* Generate control vertices for Q'' :tạo các đỉnh điều khiển cho Q'' */
	for (i = 0; i <= 1; i++) {
		Q2[i].x = (Q1[i + 1].x - Q1[i].x) * 2.0;
		Q2[i].y = (Q1[i + 1].y - Q1[i].y) * 2.0;
	}

	/* Compute Q'(u) and Q''(u)	: tính toán Q'(U) và U'' (u)*/
	Q1_u = BezierII(2, Q1, u);
	Q2_u = BezierII(1, Q2, u);

	/* Compute f(u)/f'(u) : tính toán F(u) và / f(U)' */
	numerator = (Q_u.x - P.x) * (Q1_u.x) + (Q_u.y - P.y) * (Q1_u.y);
	denominator = (Q1_u.x) * (Q1_u.x) + (Q1_u.y) * (Q1_u.y) +
		(Q_u.x - P.x) * (Q2_u.x) + (Q_u.y - P.y) * (Q2_u.y);

	/* u = u - f(u)/f'(u) */
	uPrime = u - (numerator / denominator);
	return (uPrime);
}



/*
 *  Reparameterize: xác định lại tham số
 *	Given set of points and their parameterization, try to find: tập hợp các điểm và tham số của chúng, cố gắng tìm một tham số chính xác hơn
 *   a better parameterization.
 *
 */
static double* Reparameterize(
	Point2* d,			/*  Array of digitized points : mảng các điểm số hoá	*/
	int		first, int last,		/*  Indices defining region	: chỉ số xác định vùng*/
	double* u,			/*  Current parameter values : giá trị tham số hiện tại	*/
	BezierCurve	bezCurve)	/*  Current fitted curve : đường cong thích hợp hiện tại	*/
{
	int 	nPts = last - first + 1;
	int 	i;
	double* uPrime;		/*  New parameter values: giá trị tham só mới 	*/

	uPrime = (double*)malloc(nPts * sizeof(double));
	for (i = first; i <= last; i++) {
		uPrime[i - first] = NewtonRaphsonRootFind(bezCurve, d[i], u[i -
			first]);
	}
	return (uPrime);
}






/*
 * ComputeLeftTangent, ComputeRightTangent, ComputeCenterTangent : tính toán tiếp tuyến trái, tính toán tiếp tuyến phải
 tính toán tiếp tuyến giữa
 *Approximate unit tangents at endpoints and "center" of digitized curve
 Các tiếp tuyến đơn vị gần đúng tại các điểm cuối và "tâm" của đường cong số hóa
 */
static Vector2 ComputeLeftTangent(
	Point2* d,		/*  Digitized points : điểm tham số hoá*/
	int		end)		/*  Index to "left" end of region :Chỉ mục đến "bên trái" vùng cuối */
{
	Vector2	tHat1;
	tHat1 = V2SubII(d[end + 1], d[end]);
	tHat1 = *V2Normalize(&tHat1);
	return tHat1;
}

static Vector2 ComputeRightTangent(
	Point2* d,		/*  Digitized points	: điểm số hoá	*/
	int		end)		/*  Index to "right" end of region : chỉ mục đến "bên phải"  vùng cuối */
{
	Vector2	tHat2;
	tHat2 = V2SubII(d[end - 1], d[end]);
	tHat2 = *V2Normalize(&tHat2);
	return tHat2;
}


static Vector2 ComputeCenterTangent(
	Point2* d,		/*  Digitized points	: điểm số hoá		*/
	int		center)		/*  Index to point inside region	Chỉ mục để chỉ trong khu vực */
{
	Vector2	V1, V2, tHatCenter;

	V1 = V2SubII(d[center - 1], d[center]);
	V2 = V2SubII(d[center], d[center + 1]);
	tHatCenter.x = (V1.x + V2.x) / 2.0;
	tHatCenter.y = (V1.y + V2.y) / 2.0;
	tHatCenter = *V2Normalize(&tHatCenter);
	return tHatCenter;
}



/*
 *  FitCubic :
 *  	Fit a Bezier curve to a (sub)set of digitized points :Khớp đường cong Bezier với tập hợp (phụ) các điểm được số hóa
 */
static void FitCubic(
	Point2* d,			/*  Array of digitized points */
	int		first, int last,	/* Indices of first and last pts in region */
	Vector2	tHat1, Vector2 tHat2,	/* Unit tangent vectors at endpoints */
	double	error)		/*  User-defined error squared	   */
{
	BezierCurve	bezCurve; /*Control points of fitted Bezier curve*/
	double* u;		/*  Parameter values for point :Giá trị tham số cho điểm */
	double* uPrime;	/*  Improved parameter values : Cải thiện giá trị tham số */
	double	maxError;	/*  Maximum fitting error : Lỗi sắp đăt tối đa	 */
	int		splitPoint;	/*  Point to split point set at	: Điểm để chia điểm được đặt tại */
	int		nPts;		/*  Number of points in subset : Số điểm trong tập hợp con */
	double	iterationError; /*Error below which you try iterating :Lỗi bên dưới mà bạn thử lặp lại */
	int		maxIterations = 4; /*  Max times to try iterating :Số lần tối đa để thử lặp */
	Vector2	tHatCenter;   	/* Unit tangent vector at splitPoint :Đơn vị vectơ tiếp tuyến tại điểm phân chia*/
	int		i;

	iterationError = error * error;
	nPts = last - first + 1;

	/*  Use heuristic if region only has two points in it :Sử dụng heuristic nếu khu vực chỉ có hai điểm trong đó */
	if (nPts == 2) {
		double dist = V2DistanceBetween2Points(&d[last], &d[first]) / 3.0;

		bezCurve = (Point2*)malloc(4 * sizeof(Point2));
		bezCurve[0] = d[first];
		bezCurve[3] = d[last];
		V2Add(&bezCurve[0], V2Scale(&tHat1, dist), &bezCurve[1]);
		V2Add(&bezCurve[3], V2Scale(&tHat2, dist), &bezCurve[2]);
		//	DrawBezierCurve(3, bezCurve);
		free((void*)bezCurve);
		return;
	}

	/*  Parameterize points, and attempt to fit curve :Tham số hóa các điểm và cố gắng khớp đường cong */
	u = ChordLengthParameterize(d, first, last);
	bezCurve = GenerateBezier(d, first, last, u, tHat1, tHat2);

	/*  Find max deviation of points to fitted curve :Tìm độ lệch tối đa của các điểm đối với đường cong được trang bị */
	maxError = ComputeMaxError(d, first, last, bezCurve, u, &splitPoint);
	if (maxError < error) {
		DrawBezierCurve(3, bezCurve);
		free((void*)u);
		free((void*)bezCurve);
		return;
	}


	/*  If error not too large, try some reparameterization :Nếu lỗi không quá lớn, hãy thử một số tham số lại */
	/*  and iteration  :và lặp đi lặp lại*/
	if (maxError < iterationError) {
		for (i = 0; i < maxIterations; i++) {
			uPrime = Reparameterize(d, first, last, u, bezCurve);
			bezCurve = GenerateBezier(d, first, last, uPrime, tHat1, tHat2);
			maxError = ComputeMaxError(d, first, last,
				bezCurve, uPrime, &splitPoint);
			if (maxError < error) {
				DrawBezierCurve(3, bezCurve);
				free((void*)u);
				free((void*)bezCurve);
				return;
			}
			free((void*)u);
			u = uPrime;
		}
	}

	/* Fitting failed -- split at max error point and fit recursively: sắp đặt thất bại - phân chia tại điểm lỗi tối đa và khớp đệ quy*/
	free((void*)u);
	free((void*)bezCurve);
	tHatCenter = ComputeCenterTangent(d, splitPoint);
	FitCubic(d, first, splitPoint, tHat1, tHatCenter, error);
	V2Negate(&tHatCenter);
	FitCubic(d, splitPoint, last, tHatCenter, tHat2, error);
}








/*
 *  FitCurve :
 *  	Fit a Bezier curve to a set of digitized points :Khớp đường cong Bezier với tập hợp các điểm được số hóa
 */
void FitCurve(
	Point2* d,			/*  Array of digitized points : mảng của điểm 	*/
	int		nPts,		/*  Number of digitized points	: Số điểm được số hóa*/
	double	error)		/*  User-defined error squared	: Bình phương lỗi do người dùng xác định*/
{
	Vector2	tHat1, tHat2;	/*  Unit tangent vectors at endpoints :Các vectơ tiếp tuyến đơn vị tại các điểm cuối */

	tHat1 = ComputeLeftTangent(d, 0);
	tHat2 = ComputeRightTangent(d, nPts - 1);
	FitCubic(d, 0, nPts - 1, tHat1, tHat2, error);
}


#ifdef TESTMODE

int n;


/*
 *  main:
 *	Example of how to use the curve-fitting code.  Given an array
 *   of points and a tolerance (squared error between points and
 *	fitted curve), the algorithm will generate a piecewise
 *	cubic Bezier representation that approximates the points.
 Ví dụ về cách sử dụng mã khớp đường cong. Đưa ra một mảng các điểm và dung sai (sai số bình phương giữa
các điểm và đường cong được trang bị), thuật toán sẽ tạo ra một biểu diễn Bezier khối lập phương gần đúng với các điểm.
 *	When a cubic is generated, the routine "DrawBezierCurve"
 *	is called, which outputs the Bezier curve just created
 *	(arguments are the degree and the control points, respectively).
 *	Users will have to implement this function themselves
 *   ascii output, etc.
 Khi một khối được tạo, thường trình "DrawBezierCurve" được gọi, sẽ tạo ra đường cong Bezier vừa tạo
 (các đối số lần lượt là độ và các điểm kiểm soát). Người dùng sẽ phải tự thực hiện chức năng này đầu ra ascii, v.v.
 *
 */
//void input() {
//	fstream f;
//	f.open("D:/do an/code doancoso/input.txt");
//	f >> n;
//	for (int i = 0; i < n; i++)
//		for (int j = 0; j < n; j++)
//		{
//			f >> a[i][j];
//			if (a[i][j] == 0 || i == j) a[i][j] = INFINITY;
//		}
//	f.close();
//}
int main()
{
	static Point2 d[7]=	{/*    {	  Digitized points  : tập hợp điểm*/
	{ 0.0, 0.0 },
	{ 0.0, 0.5 },
	{ 1.1, 1.4 },
	{ 2.1, 1.6 },
	{ 3.2, 1.1 },
	{ 4.0, 0.2 },
	{ 4.0, 0.0 },
}; 
	double	error = 4.0;		/*  Squared error : lỗi bình phương*/
	FitCurve(d, 7, error);		/*  Fit the Bezier curves :khớp với đường cong bezier */
}
#endif						 /* TESTMODE */
