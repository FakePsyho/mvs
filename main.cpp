// Author: Psyho
// Blog: http://psyho.gg/
// Twitter: https://twitter.com/fakepsyho

#include <thread>
#include <mutex>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <sys/time.h>

#ifdef USE_VIS
#include "CImg.h"
using namespace cimg_library;
#endif

using namespace std;

#define FOR(i,a,b)  for(int i=(a);i<(b);++i)
#define REP(i,a)    FOR(i,0,a)
#define ALL(x)      x.begin(),x.end()
#define PB          push_back
#define S           size()
#define byte        unsigned char
#define MP          make_pair
#define X           first
#define Y           second
#define VC          vector
#define PII         pair<int, int>
#define PDD         pair<double, double>
#define VI          VC<int>
#define VVI         VC<VI>
#define VPII        VC<PII>
#define VD          VC<double>
#define VVD         VC<VD>
#define VVVD        VC<VVD>
#define VS          VC<string>
#define DB(a)       cerr << #a << ": " << (a) << endl;

string NITF_DIRECTORY = "NITF";
string KML_FILE = "";
string SIMPLIFIED_DATA_DIRECTORY = "";
string OUTPUT_FILE = "";


string trim(string s) {
	int p1 = 0;
	int p2 = s.S-1;
	while (p1 < s.S && (s[p1] == ' ' || s[p1] == '\n' || s[p1] == '\r')) p1++;
	while (p2 && (s[p2] == ' ' || s[p2] == '\n' || s[p2] == '\r')) p2--;
	return p1 > p2 ? "" : s.substr(p1, p2 - p1 + 1);
}

bool startsWith(string a, string b) {
	return a.S >= b.S && a.substr(0, b.S) == b;
}

template<class A, class B> ostream& operator<<(ostream &os, pair<A,B> &p) {os << "(" << p.X << "," << p.Y << ")"; return os;}
template<class A, class B, class C> ostream& operator<<(ostream &os, tuple<A,B,C> &p) {os << "(" << get<0>(p) << "," << get<1>(p) << "," << get<2>(p) << ")"; return os;}
template<class T> ostream& operator<<(ostream &os, VC<T> &v) {os << "{"; REP(i, v.S) {if (i) os << ", "; os << v[i];} os << "}"; return os;}
template<class T> ostream& operator<<(ostream &os, set<T> &s) {VC<T> vs(ALL(s)); return os << vs;}
template<class A, class B> ostream& operator<<(ostream &os, map<A, B> &m) {VC<pair<A,B>> vs; for (auto &x : m) vs.PB(x); return os << vs;}
template<class T> string i2s(T x) {ostringstream o; if (floor(x) == x) o << (int)x; else o << x; return o.str();}
VS splt(string s, char c = ' ') {VS all; int p = 0, np; while (np = s.find(c, p), np >= 0) {all.PB(s.substr(p, np - p)); p = np + 1;} all.PB(s.substr(p)); return all;}

double getTime() {
	timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec * 1e-6;
}

VS readFile(string fn) {
	const int MAX_LINE = 1000000;
	FILE *f = fopen(fn.c_str(), "r");
	VS rv;
	while (!feof(f)) {
		static char line[MAX_LINE];
		line[0] = 0;
		fgets(line, MAX_LINE, f);
		int n = strlen(line);
		while (n && (line[n-1] == '\n' || line[n-1] == '\r')) line[--n] = 0;
		if (n == 0) continue;
		string s(line);
		s.shrink_to_fit();
		rv.PB(s);
	}
	rv.shrink_to_fit();
	return rv;
}

void saveFile(string fn, VS &data) {
	FILE *f = fopen(fn.c_str(), "w");
	for (string &s : data) {
		fputs(s.c_str(), f);
		fputs("\n", f);
	}
	fclose(f);
}

bool fileExists(string fn) {
	FILE *f = fopen(fn.c_str(), "r");
	if (f) fclose(f);
	return f != NULL;
}

string genTmpFile() {
	while (true) {
		string s = "";
		REP(i, 8) s += (char)('a' + rand() % ('z' - 'a' + 1));
		if (!fileExists(s)) return s;
	}
}

void deleteFile(string fn) {
	system(("rm " + fn).c_str());
}

VS fileList(string dir = ".") {
	string tmp = genTmpFile();
	system(("ls -1 " + dir + " > " + tmp).c_str());
	VS rv = readFile(tmp);
	REP(i, rv.S) if (rv[i] == tmp) rv.erase(rv.begin() + i);
	system(("rm " + tmp).c_str());
	return rv;
}

void system(string s) {
	system(s.c_str());
}


template <class T>
static double calcMean(VC<T> &v) {
	double rv = 0;
	REP(i, v.S) rv += v[i];
	return rv / v.S;
}

template <class T>
static double calcSTDDev(VC<T> &v) {
	double m = calcMean(v);
	double rv = 0;
	REP(i, v.S) rv += (v[i] - m) * (v[i] - m);
	return sqrt(rv / v.S);
}

double lerp(double a, double b, double t) {
	return a + (b - a) * t;
}	

int lerpi(double a, double b, double t) {
	return (int)(a + (b - a) * t);
}

int THREADS_NO = 4;

const int MAX_THREADS_NO = 4;
bool FILL_GRID = false;
double AUGMENT_DIST = 0.0;
double MERGE_AUGMENT_DIST = 0.0;
int INTEGRITY_DIST = 4;
bool EXPLOIT = false;
int MIN_OCC = 0;
int MIN2_OCC = 0;
int ERODE_DIST = 0;
double LIMIT_A = 0.2;
double LIMIT_B = 8.0;
double LIMIT_C = 100.0;
double MOVE_X = 0;
double MOVE_Y = 0;

const string RPC_PREFIX[] = {"training/rpc/rpc_Challenge1_LowRiseBuildings_CROPPED_", "testing/MasterProvisional1/rpc/rpc_MasterProvisional1_", "testing/MasterProvisional2/rpc/rpc_MasterProvisional2_", "testing/MasterProvisional3/rpc/rpc_MasterProvisional3_"};
const string KML_PREFIX[] = {"training/Challenge1.kml", "testing/MasterProvisional1/MasterProvisional1.kml", "testing/MasterProvisional2/MasterProvisional2.kml", "testing/MasterProvisional3/MasterProvisional3.kml"};

struct KML {
	double latMin = +1e9;
	double latMax = -1e9;
	double lonMin = +1e9;
	double lonMax = -1e9;
	
	static double constexpr TOLERANCE = 0.001;
	
	bool outside(double lat, double lon) {
		if (latMin > 1e8) return false;
		if (lat < latMin - TOLERANCE) return true;
		if (lat > latMax + TOLERANCE) return true;
		if (lon < lonMin - TOLERANCE) return true;
		if (lon > lonMax + TOLERANCE) return true;
		return false;
	}
};

KML readKML(string fn) {
	KML rv;
	VS vs = readFile(fn);
	FOR(i, 11, 16) {
		string s = trim(vs[i]);
		VS v = splt(s, ',');
		double lat = atof(v[0].c_str());
		double lon = atof(v[1].c_str());
		rv.latMin = min(rv.latMin, lat);
		rv.latMax = max(rv.latMax, lat);
		rv.lonMin = min(rv.lonMin, lon);
		rv.lonMax = max(rv.lonMax, lon);
	}
	return rv;
}

KML readKML(int dir) {
	return readKML(KML_PREFIX[dir]);
}

void generateFlatResult(double x1, double y1, double x2, double y2, int tx, int ty, string fn) {
	FILE *f = fopen(fn.c_str(), "w");
	REP(x, tx) REP(y, ty) {
		fprintf(f, "%.2lf %.2lf 1.0 0\n", lerp(x1, x2, 1.0 * x / (tx - 1)), lerp(y1, y2, 1.0 * y / (ty - 1)));
	}
	fclose(f);
}

double PI = 2 * acos(0);
PDD convertCoords2UTM(double lat, double lon) {
	// based on http://stackoverflow.com/questions/176137/java-convert-lat-lon-to-utm
	
	int zone = (int) floor(lon/6+31);
	int letterIndex = (int)(lat + 80) / 8;
	if (letterIndex < 0) letterIndex = 0;
	if (letterIndex >= 20) letterIndex = 19;
	// char letter = "CDEFGHJKLMNPQRSTUVWX".charAt(letterIndex) + "";
	
	double DEG = PI / 180;
	
	double latR = lat * DEG;
	double lonR = lon * DEG;
	
	double easting = 0.5 * log((1 + cos(latR) * sin(lonR - (6 * zone - 183) * DEG)) / (1 - cos(latR) * sin(lonR - (6 * zone - 183) * DEG))) * 0.9996 * 6399593.62 / pow((1 + pow(0.0820944379, 2) * pow(cos(latR), 2)), 0.5) * (1 + pow(0.0820944379, 2) / 2 * pow((0.5 * log((1 + cos(latR) * sin(lonR - (6 * zone - 183) * DEG)) / (1 - cos(latR) * sin(lonR - (6 * zone - 183) * DEG)))), 2) * pow(cos(latR), 2) / 3) + 500000;
	double northing = (atan(tan(latR) / cos((lonR - (6 * zone - 183) * DEG))) - latR) * 0.9996 * 6399593.625 / sqrt(1 + 0.006739496742 * pow(cos(latR), 2)) * (1 + 0.006739496742 / 2 * pow(0.5 * log((1 + cos(latR) * sin((lonR - (6 * zone - 183) * DEG))) / (1 - cos(latR) * sin((lonR - (6 * zone - 183) * DEG)))), 2) * pow(cos(latR), 2)) + 0.9996 * 6399593.625 * (latR - 0.005054622556 * (latR + sin(2 * latR) / 2) + 4.258201531e-05 * (3 * (latR + sin(2 * latR) / 2) + sin(2 * latR) * pow(cos(latR), 2)) / 4 - 1.674057895e-07 * (5 * (3 * (latR + sin(2 * latR) / 2) + sin(2 * latR) * pow(cos(latR), 2)) / 4 + sin(2 * latR) * pow(cos(latR), 2) * pow(cos(latR), 2)) / 3);
	if (lat < -8)
		northing = northing + 10000000;
	return MP(easting, northing);
}

map<string, string> readMeta(string fn) {
	map<string, string> m;
	VS v = readFile(fn);
	for (string &s : v) {
		VS vs = splt(s, '\t');
		assert(vs.S == 1 || vs.S == 2);
		if (vs.S == 1) vs.PB("");
		m[vs[0]] = vs[1];
	}
	return m;
}

void convertRPC2XML(string fin, string fout) {
	map<string, string> m = readMeta(fin);
	VS v;
	v.PB("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>");
	v.PB("<isd>");
	v.PB("<RPB>");
	v.PB("<SATID>WV03</SATID>");
	v.PB("<BANDID>MS1</BANDID>");
	v.PB("<SPECID>RPC00B</SPECID>");
	v.PB("<IMAGE>");
	VS fields = {"ERRBIAS","RPC00B_ERR_BIAS","ERRRAND","RPC00B_ERR_RAND","LINEOFFSET","RPC00B_LINE_OFF","SAMPOFFSET","RPC00B_SAMP_OFF","LATOFFSET","RPC00B_LAT_OFF","LONGOFFSET","RPC00B_LONG_OFF","HEIGHTOFFSET","RPC00B_HEIGHT_OFF","LINESCALE","RPC00B_LINE_SCALE","SAMPSCALE","RPC00B_SAMP_SCALE","LATSCALE","RPC00B_LAT_SCALE","LONGSCALE","RPC00B_LONG_SCALE","HEIGHTSCALE","RPC00B_HEIGHT_SCALE"};
	VS arrays = {"LINENUMCOEF", "RPC00B_LINE_NUM_COEFF", "LINEDENCOEF", "RPC00B_LINE_DEN_COEFF", "SAMPNUMCOEF", "RPC00B_SAMP_NUM_COEFF", "SAMPDENCOEF", "RPC00B_SAMP_DEN_COEFF"};
	for (int i = 0; i < fields.S; i += 2)
		v.PB("<" + fields[i] + ">" + m[fields[i+1]] + "</" + fields[i] + ">");
	for (int i = 0; i < arrays.S; i += 2) {
		string s = "";
		s += "<" + arrays[i] + "List><" + arrays[i] + ">";
		REP(j, 20) {
			if (j) s += " ";
			s += m[arrays[i+1]+"_"+i2s(j+1)];
		}
		s += "</" + arrays[i] + "></" + arrays[i] + "List>";
		v.PB(s);
	}
	
	v.PB("</IMAGE>");
	v.PB("</RPB>");
	v.PB("</isd>");
	saveFile(fout, v);
}

PDD readRPCOffset(string dir, string file) {
	string fn = dir + file + ".txt";
	if (!fileExists(fn)) return MP(-1, -1);
	VS v = splt(readFile(fn)[0], ',');
	return MP(atof(v[94].c_str()), atof(v[95].c_str()));
}

PDD readRPCOffset(int dir, string file) {
	return readRPCOffset(RPC_PREFIX[dir], file);
}

void convertTruth2UTM(string fin, string fout, int pruneDir = -1, string kml_file = "") {
	VS v1 = readFile(fin);
	VS v2;
	
	KML kml;
	if (pruneDir != -1 && kml_file.S == 0) kml = readKML(pruneDir);
	if (kml_file.S) kml = readKML(kml_file);
	for (string s : v1) {
		VS v = splt(s, ',');
		double lat = atof(v[0].c_str());
		double lon = atof(v[1].c_str());
		lat += MOVE_X;
		lon += MOVE_Y;
		if (kml.outside(lat, lon)) continue;
		PDD p = convertCoords2UTM(lon, lat);
		char str[200];
		if (v.S <= 3) v.PB("0");
		sprintf(str, "%.3lf %.3lf %s %s", p.X, p.Y, v[2].c_str(), v[3].c_str());
		v2.PB(str);
	}
	saveFile(fout, v2);
}

struct Pos {
	double x;
	double y;
	double h;
	Pos() {
		x = y = h = 0;
	}
	
	Pos(double _x, double _y, double _h) {
		x = _x;
		y = _y;
		h = _h;
	}
	
	bool empty() {
		return x == 0 && y == 0 && h == 0.0;
	}
	
	friend ostream& operator<< (ostream &os, Pos &p) {
		os << "(" << p.x << "," << p.y << "," << p.h << ")"; return os;
	}
};

double medianArray[3000 * 3000];
int medianArrayNo = 0;
double ma[MAX_THREADS_NO][3000 * 3000];
int mno[MAX_THREADS_NO];


struct Grid {
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double block;
	VVD v;
	VVI vt;
	VVVD va;
	int xs = 0;
	int ys = 0;
	bool empty = true;
	
	void initGrid(KML kml, double _block) {
		block = _block;
		xmin = +1e9;
		xmax = -1e9;
		ymin = +1e9;
		ymax = -1e9;
		REP(i, 2) REP(j, 2) {
			double lat = i == 0 ? kml.latMin : kml.latMax;
			double lon = j == 0 ? kml.lonMin : kml.lonMax;
			PDD p1 = convertCoords2UTM(lon - KML::TOLERANCE, lat - KML::TOLERANCE);
			PDD p2 = convertCoords2UTM(lon + KML::TOLERANCE, lat + KML::TOLERANCE);
			xmin = min(xmin, min(p1.X, p2.X));
			xmax = max(xmax, max(p1.X, p2.X));
			ymin = min(ymin, min(p1.Y, p2.Y));
			ymax = max(ymax, max(p1.Y, p2.Y));
		}
		xs = (int)floor((xmax - xmin) / block) + 1;
		ys = (int)floor((ymax - ymin) / block) + 1;
		v = VVD(xs, VD(ys, 0.0));
		vt = VVI(xs, VI(ys, 0));
		va = VVVD(xs, VVD(ys, VD()));
	}
	
	void initGrid(Grid &grid) {
		xmin = grid.xmin;
		xmax = grid.xmax;
		ymin = grid.ymin;
		ymax = grid.ymax;
		block = grid.block;
		xs = grid.xs;
		ys = grid.ys;
		v = VVD(xs, VD(ys, 0.0));
		vt = VVI(xs, VI(ys, 0));
		va = VVVD(xs, VVD(ys, VD()));
	}
	
	void clearGrid() {
		REP(x, xs) REP(y, ys) v[x][y] = 0;
		REP(x, xs) REP(y, ys) vt[x][y] = 0;
		// REP(x, xs) REP(y, ys) va[x][y].clear();
	}
	
	void clearGridFull() {
		REP(x, xs) REP(y, ys) v[x][y] = 0;
		REP(x, xs) REP(y, ys) vt[x][y] = 0;
		REP(x, xs) REP(y, ys) va[x][y].clear();
	}
	
	PII getPos(double x, double y) {
		return MP((int)floor((x - xmin) / block), (int)floor((y - ymin) / block));
	}
	
	void addGrid(VC<Pos> &vp, Pos offset = Pos()) {
		assert(xs > 0 && ys > 0);
		for (Pos &pos : vp) {
			PII p = getPos(pos.x + offset.x, pos.y + offset.y);
			if (p.X < 0 || p.X >= xs || p.Y < 0 || p.Y >= ys) continue;
			v[p.X][p.Y] += pos.h + offset.h;
			vt[p.X][p.Y]++;
			va[p.X][p.Y].PB(pos.h + offset.h);
		}
		empty = false;
	}
	
	void addGrid(Grid &g) {
		assert(xs == g.xs && ys == g.ys);
		REP(x, xs) REP(y, ys) {
			v[x][y] += g.v[x][y];
			vt[x][y] += g.vt[x][y];
			for (double &d : g.va[x][y]) va[x][y].PB(d);
		}
		empty = false;
	}
	
	void removeGrid(Grid &g) {
		assert(xs == g.xs && ys == g.ys);
		REP(x, xs) REP(y, ys) {
			v[x][y] -= g.v[x][y];
			vt[x][y] -= g.vt[x][y];
			for (double &d : g.va[x][y]) {
				auto it = find(ALL(va[x][y]), d);
				assert(it != va[x][y].end());
				va[x][y].erase(it);
			}
		}
		empty = false;
	}
	
	void simplifyGrid() {
		REP(x, xs) REP(y, ys) if (vt[x][y] > 1) {
			v[x][y] /= vt[x][y];
			vt[x][y] = 1;
		}
		if (FILL_GRID) {
			FOR(x, 1, xs-1) FOR(y, 1, ys-1) if (vt[x][y] == 0 && vt[x-1][y] && vt[x+1][y] && vt[x][y-1] && vt[x][y+1]) {
				v[x][y] = (v[x-1][y] + v[x+1][y] + v[x][y-1] + v[x][y+1]) / 4;
				vt[x][y] = 1;
				va[x][y].PB(v[x][y]);
			}
		}
		
		if (ERODE_DIST) {
			VVD tv = v;
			REP(x, xs) REP(y, ys) {
				if (!vt[x][y]) continue;
				FOR(dx, -ERODE_DIST, ERODE_DIST + 1) FOR(dy, -ERODE_DIST, ERODE_DIST+1) {
					if (abs(dx) + abs(dy) > ERODE_DIST) continue;
					int nx = x + dx;
					int ny = y + dy;
					if (nx < 0 || nx >= xs || ny < 0 || ny >= ys) continue;
					tv[x][y] = min(tv[x][y], v[nx][ny]);
				}
			}
			v = tv;
		}
		// REP(x, xs) REP(y, ys) if (vt[x][y]) va[x][y].PB(v[x][y]);
	}
	
	void addGridTmp(VC<Pos> &vp, Pos offset = Pos()) {
		assert(xs > 0 && ys > 0);
		for (Pos &pos : vp) {
			PII p = getPos(pos.x + offset.x, pos.y + offset.y);
			if (p.X < 0 || p.X >= xs || p.Y < 0 || p.Y >= ys) continue;
			v[p.X][p.Y] += pos.h + offset.h;
			vt[p.X][p.Y]++;
		}
		empty = false;
	}
	
	void updateGridTmp(double diff) {
		REP(x, xs) REP(y, ys) v[x][y] += diff * vt[x][y];
	}
	
	double calcDiff(Grid &g) {
		double sum = 0;
		int cnt = 0;
		REP(x, xs) REP(y, ys) if (vt[x][y] && g.vt[x][y]) {
			double d = v[x][y] / vt[x][y] - g.v[x][y] / g.vt[x][y];
			sum += d;
			cnt++;
		}
		return sum / cnt;
	}
	
	double calcXDiff(Grid &g, int a = 0) {
		mno[a] = 0;
		REP(x, xs) REP(y, ys) if (vt[x][y] && g.vt[x][y]) {
			double d = v[x][y] / vt[x][y] - g.v[x][y] / g.vt[x][y];
			ma[a][mno[a]++] = d;
		}
		sort(ma[a], ma[a] + mno[a]);
		int cnt = 0;
		double sum = 0;
		FOR(i, (int)(mno[a] * 0.2), (int)(mno[a] * 0.8)) {
			sum += ma[a][i];
			cnt++;
		}
		return sum / cnt;
	}
	
	double calcMedian(Grid &g, int a = 0) {
		mno[a] = 0;
		REP(x, xs) REP(y, ys) if (vt[x][y] && g.vt[x][y]) {
			double d = v[x][y] / vt[x][y] - g.v[x][y] / g.vt[x][y];
			ma[a][mno[a]++] = d;
		}
		nth_element(ma[a], ma[a] + mno[a] / 2, ma[a] + mno[a]);
		return ma[a][mno[a]/2];
	}
	
	double calcMSE(Grid &g) {
		double sum = 0;
		int cnt = 0;
		REP(x, xs) REP(y, ys) if (vt[x][y] && g.vt[x][y]) {
			double d = v[x][y] / vt[x][y] - g.v[x][y] / g.vt[x][y];
			sum += d * d;
			cnt++;
		}
		return sqrt(sum / cnt);
	}
	
	double calcMAE(Grid &g) {
		double sum = 0;
		int cnt = 0;
		REP(x, xs) REP(y, ys) if (vt[x][y] && g.vt[x][y]) {
			double d = v[x][y] / vt[x][y] - g.v[x][y] / g.vt[x][y];
			sum += abs(d);
			cnt++;
		}
		return sum / cnt;
	}
	
	Pos findOffset(VC<Pos> &vp, Pos initialOffset = Pos(0, 0, 0.0), bool optimize = false) {
		Pos best(0.0, 0.0, 0.0);
		double bv = 1e30;
		
		Grid g[THREADS_NO];
		REP(i, THREADS_NO) g[i].initGrid(*this);
		
		double step = 0.3 * 32;
		int level = 0;
		
		if (optimize) {
			step = 0.3 * 2;
			best = initialOffset;
		}
		
		while (true) {
			Pos center = best;
			VPII vv;
			int vpos = 0;
			FOR(dx, -1, 2) FOR(dy, -1, 2) 
				if (dx || dy || level == 0)
					vv.PB(MP(dx, dy));
			mutex mut;
			thread *threads = new thread[THREADS_NO];
			REP(t, THREADS_NO) {
				threads[t] = thread([&](int i) {
					while (true) {
						mut.lock();
						int pos = vpos++;
						mut.unlock();
						if (pos >= vv.S) break;
						int dx = vv[pos].X;
						int dy = vv[pos].Y;
						
						g[i].clearGrid();
						Pos offset = Pos(center.x + dx * step, center.y + dy * step, 0);
						g[i].addGridTmp(vp, offset);
						double diff = optimize ? calcXDiff(g[i], i) : calcMedian(g[i], i);
						offset.h = diff;
						g[i].updateGridTmp(diff);
						
						double av = calcMSE(g[i]);
						mut.lock();
						if (av < bv) {
							best = offset;
							bv = av;
						}
						mut.unlock();
					}
				}, t);
			}
			REP(i, THREADS_NO) threads[i].join();
			delete[] threads;
			if (step < 0.1) break;
			level++;
			step /= 2;
		}
		cout << "Best Offset: " << best.x << ' ' << best.y << ' ' << best.h << ' ' << bv << endl;
		return best;
	}
	
	Pos mergeGrid(VC<Pos> &vp, Pos originalOffset = Pos(0, 0, 0.0), bool final = false, bool optimize = false) {
		if (final) {
			Grid g;
			g.initGrid(*this);
			g.addGrid(vp, originalOffset);
			g.simplifyGrid();
			addGrid(g);
			return originalOffset;
		}
		
		if (empty) {
			Pos offset = Pos();
			addGrid(vp, offset);
			simplifyGrid();
			return offset;
		} else {
			Pos offset = findOffset(vp, originalOffset, optimize);
			Grid g;
			g.initGrid(*this);
			g.addGrid(vp, offset);
			g.simplifyGrid();
			addGrid(g);
			return offset;
		}
	}
	
	void removeGrid(VC<Pos> &vp, Pos offset) {
		Grid g;
		g.initGrid(*this);
		g.addGrid(vp, offset);
		g.simplifyGrid();
		removeGrid(g);
	}
	
	void save(string fout) {
		FILE *f = fopen(fout.c_str(), "w");
		VVD nv(xs, VD(ys, 0));
		REP(x, xs) REP(y, ys) if (vt[x][y]) {
			sort(ALL(va[x][y]));
			if (va[x][y].S < 5) {
				nv[x][y] = va[x][y].S % 2 == 1 ? va[x][y][va[x][y].S/2] : (va[x][y][va[x][y].S/2-1] + va[x][y][va[x][y].S/2]) / 2;
			} else {
				int p0 = 0;
				int p1 = va[x][y].S;
				while (p1 - p0 > 5) {
					p0++;
					p1--;
				}
				double sum = 0;
				int cnt = 0;
				FOR(i, p0, p1) {
					sum += va[x][y][i];
					cnt++;
				}
				nv[x][y] = sum / cnt;
			}
		}
		int DD = INTEGRITY_DIST/2;
		REP(x, xs) REP(y, ys) if (vt[x][y]) {
			bool useFilter = false;
			
			bool ok = false;
			if (EXPLOIT) {
				double mn = +1e9;
				double mx = -1e9;
				int tot = 0;
				bool bad = false;
				FOR(dx, -DD, DD+1) FOR(dy, -DD, DD+1) {
					int bigd = max(abs(dx), abs(dy));
					int smalld = min(abs(dx), abs(dy));
					if (bigd*2+smalld>INTEGRITY_DIST) continue;
					int nx = x + dx;
					int ny = y + dy;
					if (nx < 0 || nx >= xs || ny < 0 || ny >= ys)
						continue;
					if (vt[nx][ny]) {
						mx = max(mx, nv[nx][ny]);
						mn = min(mn, nv[nx][ny]);
						tot += vt[nx][ny];
					} else {
						bad = true;
					}
				}
				if (vt[x][y] == 0) continue;
				double rmse = calcSTDDev(va[x][y]);
				ok = false;
				int cno = 0;
				FOR(dx, -1, 2) FOR(dy, -1, 2) {
					if (dx && dy) continue;
					int nx = x + dx;
					int ny = y + dy;
					if (nx < 0 || nx >= xs || ny < 0 || ny >= ys)
						continue;
					cno += va[x][y].S;
				}
				if (abs(mn - mx) < LIMIT_A && tot >= MIN_OCC && !bad && cno > 5) ok = true;			
				if (x%3==0&&y%3==0&&!ok&&abs(mn-mx)<LIMIT_B&&tot>=MIN2_OCC && cno > 5) ok = true, useFilter = true;
			} else {
				ok = (x+y)%2;
			}
			
			if (!ok) continue;
			
			double height = 0;
			if (useFilter) {
				medianArrayNo = 0;
				FOR(dx, -1, 2) FOR(dy, -1, 2) {
					int nx = x + dx;
					int ny = y + dy;
					if (nx < 0 || nx >= xs || ny < 0 || ny >= ys || vt[nx][ny] == 0) continue;
					REP(i, vt[nx][ny]) medianArray[medianArrayNo++] = va[nx][ny][i];
				}
				sort(medianArray, medianArray + medianArrayNo);
				height = medianArrayNo % 2 ? medianArray[medianArrayNo/2] : (medianArray[medianArrayNo/2-1] + medianArray[medianArrayNo/2]) / 2;
			} else {
				height = nv[x][y];
			}
			
			fprintf(f, "%.3lf %.3lf %.4lf 0\n", xmin + (x + 0.5) * block, ymin + (y + 0.5) * block, height);
		}
		fclose(f);
	}
};

VC<Pos> readTruth(string fin) {
	VC<Pos> rv;
	VS vs = readFile(fin);
	for (string &s : vs) {
		VS va = splt(s, ' ');
		assert(va.S >= 3);
		rv.PB(Pos(atof(va[0].c_str()), atof(va[1].c_str()), atof(va[2].c_str())));
	}
	return rv;
}

void augmentTruth(VC<Pos> &v, double d) {
	int sz = v.S;
	REP(i, sz) {
		v.PB(Pos(v[i].x + d, v[i].y + 0, v[i].h));
		v.PB(Pos(v[i].x - d, v[i].y + 0, v[i].h));
		v.PB(Pos(v[i].x + 0, v[i].y + d, v[i].h));
		v.PB(Pos(v[i].x + 0, v[i].y - d, v[i].h));
	}
}

void sortNITFFiles(VS &vs) {
	VS months = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
	map<string, int> monthMap;
	REP(i, months.S) monthMap[months[i]] = i;
	sort(ALL(vs), [&](const string &a, const string &b) -> bool {
		int day1 = atoi(a.substr(0, 2).c_str());
		int day2 = atoi(b.substr(0, 2).c_str());
		int month1 = monthMap[a.substr(2, 3)];
		int month2 = monthMap[b.substr(2, 3)];
		int year1 = atoi(a.substr(5, 2).c_str());
		int year2 = atoi(b.substr(5, 2).c_str());
		if (year1 != year2) return year1 < year2;
		if (month1 != month2) return month1 < month2;
		if (day1 != day2) return day1 < day2;
		return a < b;
	});
}

pair<string, string> getNITFPair(int id) {
	VS files = fileList(NITF_DIRECTORY);
	for (string &s : files) s = s.substr(0, s.S - 4);
	sortNITFFiles(files);
	
	int diff = 1;
	int a = -1;
	while (true) {
		a++;
		int b = a + diff;
		if (b >= files.S) {
			a = -1;
			diff++;
			continue;
		}
		id--;
		if (id < 0) return MP(files[a], files[b]);
	}
}



void visualizeTruth(string fin, string fout, double scale = 1) {
#ifdef USE_VIS
	VC<tuple<double, double, double>> v;
	VS vs = readFile(fin);
	DB(vs.S);
	for (string &s : vs) {
		VS va = splt(s, ' ');
		assert(va.S >= 3);
		v.PB(make_tuple(atof(va[0].c_str()), atof(va[1].c_str()), atof(va[2].c_str())));
	}
	DB(v.S);
	double maxx = -1e9, minx = 1e9;
	double maxy = -1e9, miny = 1e9;
	VD vh;
	for (auto &t : v) {
		maxx = max(maxx, get<0>(t));
		minx = min(minx, get<0>(t));
		maxy = max(maxy, get<1>(t));
		miny = min(miny, get<1>(t));
		if (get<2>(t) != -9999)
			vh.PB(get<2>(t));
	}
	sort(ALL(vh));
	
	DB(minx); DB(maxx);
	DB(miny); DB(maxy);
	int width = (maxx - minx) / 0.30 + 1;
	int height = (maxy - miny) / 0.30 + 1;
	minx -= 0.15;
	miny -= 0.15;
	
	double ground = vh[(int)(vh.S * 0.1)];
	
	VVD sum(width, VD(height));
	VVD tot(width, VD(height));
	for (auto &t : v) {
		if (get<2>(t) == -9999) continue;
		int x = min(width-1, (int)((get<0>(t) - minx) / (maxx - minx) * width));
		int y = min(height-1, (int)((get<1>(t) - miny) / (maxy - miny) * height));
		assert(x >= 0 && x < width);
		assert(y >= 0 && y < height);
		double h = get<2>(t) - ground;
		sum[x][y] += h;
		tot[x][y] += 1;
	}
	
	VVI colors = {{107, 47, 107}, {238, 100, 238}, {64, 64, 255}, {0, 255, 255}, {64, 255, 64}, {255, 200, 0}, {255, 64, 64}};
	VD heights = {-30, -10, 0, 5, 15, 30, 50};
	
	CImg<byte> img(width, height, 1, 3);
	img = (byte)0;
	REP(x, width) REP(y, height) {
		if (tot[x][y] == 0) continue;
		double h = sum[x][y] / tot[x][y];
		h *= scale;
		if (h < heights[0]) REP(i, 3) img(x, y, 0, i) = colors[0][i];
		if (h > heights.back()) REP(i, 3) img(x, y, 0, i) = colors.back()[i];
		FOR(j, 1, heights.S) if (h >= heights[j-1] && h <= heights[j]) REP(i, 3) img(x, y, 0, i) = lerpi(colors[j-1][i], colors[j][i], (h - heights[j-1]) / (heights[j] - heights[j-1]));
	}
	img.save(fout.c_str());
#else
	cout << "[Error] USE_VIS not defined" << endl;
#endif
}

int main(int argc, char **argv) {
	VS arguments;
	REP(i, argc) arguments.PB(argv[i]);
	DB(arguments);
	
	string visInput = "", visOutput = "";
	string convertInput = "", convertOutput = "";
	string prepareMode = "";
	string input = "";
	string output = "";
	string options = "";
	string parameters = "";
	VS mergeFiles;
	VI idList;
	
	double visScale = 1.0;
	VD mergeOffsets;
	
	bool optimize = false;
	
	bool runMode = false;
	bool lasMode = false;
	bool utmMode = false;
	bool deleteTmpFiles = false;
	bool scoreMode = false;
	bool visualizeScoreMode = false;
	bool randomMatchMode = false;
	bool mergeMode = false;
	bool snapMode = false;
	bool finalMode = false;
	
	bool projectionEnabled = false;
	
	int dir = -1;
	string fileA = "";
	string fileB = "";
	string outputPrefix = "";
	
	FOR(i, 1, argc) {
		string cmd = argv[i];
		if (cmd == "-v") {
			visInput = argv[++i];
			if (i+1 >= argc || argv[i+1][0] == '-')
				visOutput = splt(visInput, '.')[0] + ".bmp";
			else 
				visOutput = argv[++i];
		} else if (cmd == "-off") {
			while (i+1 < argc && (argv[i+1][0] != '-' || isdigit(argv[i+1][1])))
				mergeOffsets.PB(atof(argv[++i]));
		} else if (cmd == "-t") {
			THREADS_NO = atoi(argv[++i]);
		} else if (cmd == "-opt") {
			optimize = true;
		} else if (cmd == "-sv") {
			visualizeScoreMode = true;
			visInput = argv[++i];
		} else if (cmd == "-fill") {
			FILL_GRID = true;
		} else if (cmd == "-exp") {
			EXPLOIT = true;
		} else if (cmd == "-aug") {
			AUGMENT_DIST = atof(argv[++i]);
		} else if (cmd == "-maug") {
			MERGE_AUGMENT_DIST = atof(argv[++i]);
		} else if (cmd == "-move") {
			MOVE_X = atof(argv[++i]);
			MOVE_Y = atof(argv[++i]);
		} else if (cmd == "-occ") {
			MIN_OCC = atoi(argv[++i]);
		} else if (cmd == "-occ2") {
			MIN2_OCC = atoi(argv[++i]);
		} else if (cmd == "-check") {
			INTEGRITY_DIST = atoi(argv[++i]);
		} else if (cmd == "-erode") {
			ERODE_DIST = atoi(argv[++i]);
		} else if (cmd == "-lima") {
			LIMIT_A = atof(argv[++i]);
		} else if (cmd == "-limb") {
			LIMIT_B = atof(argv[++i]);
		} else if (cmd == "-limc") {
			LIMIT_C = atof(argv[++i]);
		} else if (cmd == "-seed") {
			srand(atoi(argv[++i]));
		} else if (cmd == "-project") {
			projectionEnabled = true;
		} else if (cmd == "-vscale") {
			visScale = atof(argv[++i]);
		} else if (cmd == "-cv") {
			utmMode = true;
			convertInput = argv[++i];
			convertOutput = argv[++i];
		} else if (cmd == "-lv") {
			lasMode = true;
			utmMode = true;
			convertInput = argv[++i];
			convertOutput = argv[++i];
		} else if (cmd == "-p") {
			prepareMode = argv[++i];
		} else if (cmd == "-f") {
			finalMode = true;
			dir = atoi(argv[++i]);
			outputPrefix = argv[++i];
			while (i+1 < argc && argv[i+1][0] != '-') {
				i++;
				if (argv[i][0] == '+') {
					options = &(argv[i][1]);
				} else if (argv[i][0] == '=') {
					parameters = &(argv[i][1]);
				} else {
					idList.PB(atoi(argv[i]));
				}
			}
		} else if (cmd == "-r") {
			runMode = true;
			dir = atoi(argv[++i]);
			outputPrefix = argv[++i];
			fileA = argv[++i];
			fileB = argv[++i];
			if (i+1 < argc && argv[i+1][0] == '+') 
				options = &(argv[++i][1]);
			if (i+1 < argc && argv[i+1][0] == '=') 
				parameters = &(argv[++i][1]);
		} else if (cmd == "-nitfdir") {
			NITF_DIRECTORY = argv[++i];
		} else if (cmd == "-kml") {
			KML_FILE = argv[++i];
		} else if (cmd == "-data") {
			SIMPLIFIED_DATA_DIRECTORY = argv[++i];
			while (SIMPLIFIED_DATA_DIRECTORY.S && SIMPLIFIED_DATA_DIRECTORY[SIMPLIFIED_DATA_DIRECTORY.S-1] == '/') SIMPLIFIED_DATA_DIRECTORY = SIMPLIFIED_DATA_DIRECTORY.substr(0, SIMPLIFIED_DATA_DIRECTORY.S - 1);
		} else if (cmd == "-output") {
			OUTPUT_FILE = argv[++i];
		} else if (cmd == "-d") {
			deleteTmpFiles = true;
		} else if (cmd == "-m") {
			mergeMode = true;
			dir = atoi(argv[++i]);
			while (i+2 < argc && argv[i+2][0] != '-')
				mergeFiles.PB(argv[++i]);
			output = argv[++i];
		} else if (cmd == "-s") {
			scoreMode = true;
			input = argv[++i];
		} else if (cmd == "-snap") {
			snapMode = true;
			dir = atoi(argv[++i]);
			input = argv[++i];
			output = argv[++i];
		} else if (cmd == "-x") {
			randomMatchMode = true;
			output = argv[++i];
			if (i+1 < argc && argv[i+1][0] == '+') 
				options = &(argv[++i][1]);
			if (i+1 < argc && argv[i+1][0] == '=') 
				parameters = &(argv[++i][1]);
		} else {
			cout << "[Error] Unknown Command: " << cmd << endl;
			exit(1);
		}
	}
	
	string RPCDir = SIMPLIFIED_DATA_DIRECTORY.S ? SIMPLIFIED_DATA_DIRECTORY + "/rpc/rpc_" + splt(SIMPLIFIED_DATA_DIRECTORY, '/').back() + "_" : RPC_PREFIX[dir];
	
	if (finalMode) {
		VS resultFiles;
		for (int id : idList) {
			auto files = getNITFPair(id);
			if (readRPCOffset(RPCDir, files.X).X == -1) continue;
			if (readRPCOffset(RPCDir, files.Y).X == -1) continue;
			char idstr[10];
			sprintf(idstr, "%04d", id);
			string rf = outputPrefix + string(idstr);
			resultFiles.PB(rf);
			system("./x -nitfdir " + NITF_DIRECTORY + " -kml " + KML_FILE + " -data " + SIMPLIFIED_DATA_DIRECTORY + " -output " + OUTPUT_FILE + " -r " + i2s(dir) + " " + rf + " " + files.X + " " + files.Y + " \"+" + options + "\" \"=" + parameters + "\"" + (projectionEnabled ? string(" -project") : string()));
			#ifdef USE_VIS
			system("./x -v results/" + rf + ".txt " + rf + ".bmp");
			#endif
		}
	}
	
	if (visualizeScoreMode) {
		VC<Pos> pred = readTruth(visInput);
		VC<Pos> truth = readTruth(visInput);
	}
	
	if (snapMode) {
		KML kml = KML_FILE == "" ? readKML(dir) : readKML(KML_FILE);
		VC<Pos> v = readTruth(input);
		Grid grid;
		grid.initGrid(kml, 0.30);
		grid.addGrid(v);
		grid.save(output);
	}
	
	if (mergeMode) {
		KML kml = KML_FILE == "" ? readKML(dir) : readKML(KML_FILE);
		Grid grid;
		grid.initGrid(kml, 0.30);
		VC<Pos> offsets;
		VC<VC<Pos>> truths;
		
		if (mergeOffsets.S == 0) {
			for (string &input : mergeFiles) {
				VC<Pos> v = readTruth(input);
				if (MERGE_AUGMENT_DIST) augmentTruth(v, MERGE_AUGMENT_DIST);
				truths.PB(v);
				Pos offset = grid.mergeGrid(v);
				offsets.PB(offset);
			}
			REP(i, mergeFiles.S) {
				string &input = mergeFiles[i];
				grid.removeGrid(truths[i], offsets[i]);
				offsets[i] = grid.mergeGrid(truths[i], offsets[i], false, optimize);
			}
			
			cout << "Offsets:";
			REP(i, offsets.S) cout << " " << offsets[i].x << ' ' << offsets[i].y << ' ' << offsets[i].h;
			cout << endl;
			
			if (MERGE_AUGMENT_DIST) {
				grid.clearGridFull();
				REP(i, mergeFiles.S) {
					VC<Pos> v = readTruth(mergeFiles[i]);
					if (AUGMENT_DIST) augmentTruth(v, AUGMENT_DIST);
					grid.mergeGrid(v, offsets[i], true);
				}
			}
			
		} else {
			int pos = 0;
			assert(mergeOffsets.S == mergeFiles.S * 3);
			for (string &input : mergeFiles) {
				VC<Pos> v = readTruth(input);
				if (AUGMENT_DIST) augmentTruth(v, AUGMENT_DIST);
				grid.mergeGrid(v, Pos(mergeOffsets[pos+0], mergeOffsets[pos+1], mergeOffsets[pos+2]), true);
				pos += 3;
			}
		}
		grid.save(output);
	}
	
	if (visInput.S) {
		cout << "Visualizing Truth FIle" << endl;
		visualizeTruth(visInput, visOutput, visScale);
	}
	
	if (lasMode) {
		cout << "Converting Las to XYZ" << endl;
		system(("./las2txt.exe -i " + convertInput + " -parse xyz").c_str());
		if (deleteTmpFiles) deleteFile(convertInput);
		convertInput = convertInput.substr(0, convertInput.S - 4) + ".txt";
	}
	
	if (utmMode) {
		cout << "Converting Truth File to ULM format" << endl;
		convertTruth2UTM(convertInput, convertOutput);
		if (deleteTmpFiles) deleteFile(convertInput);
	}
	
	if (runMode) {
		PDD p1 = readRPCOffset(RPCDir, fileA);
		PDD p2 = readRPCOffset(RPCDir, fileB);
		assert(p1.X != -1 && p2.X != -1);
		deleteFile("stereo.default");
		if (!projectionEnabled) {
			options += " --left-image-crop-win " + i2s((int)p1.X) + " " + i2s((int)p1.Y) + (dir == 0 ? " 3001 3001" : " 2001 2001");
			options += " --right-image-crop-win " + i2s((int)p2.X) + " " + i2s((int)p2.Y) + (dir == 0 ? " 3001 3001" : " 2001 2001");
		}
		
		if (projectionEnabled) {
			KML kml = KML_FILE == "" ? readKML(dir) : readKML(KML_FILE);
			system("mapproject -t rpc dem.tif " + NITF_DIRECTORY + "/" + fileA + ".NTF tmp1.tif --t_projwin " + i2s(kml.latMin - KML::TOLERANCE) + " " + i2s(kml.lonMin - KML::TOLERANCE) + " " + i2s(kml.latMax + KML::TOLERANCE) + " " + i2s(kml.lonMax + KML::TOLERANCE));
			system("mapproject -t rpc dem.tif " + NITF_DIRECTORY + "/" + fileB + ".NTF tmp2.tif --t_projwin " + i2s(kml.latMin - KML::TOLERANCE) + " " + i2s(kml.lonMin - KML::TOLERANCE) + " " + i2s(kml.latMax + KML::TOLERANCE) + " " + i2s(kml.lonMax + KML::TOLERANCE));
			system("stereo -t rpcmaprpc --alignment-method none tmp1.tif tmp2.tif " + NITF_DIRECTORY + "/" + fileA + ".NTF " + NITF_DIRECTORY + "/" + fileB + ".NTF results/" + outputPrefix + " dem.tif " + options);
		} else {
			system("stereo " + NITF_DIRECTORY + "/" + fileA + ".NTF " + NITF_DIRECTORY + "/" + fileB + ".NTF results/" + outputPrefix + " " + options);
		}
		system("point2las results/" + outputPrefix + "-PC.tif");
		system("las2txt -i results/" + outputPrefix + "-PC.las --parse xyz -o results/" + outputPrefix + ".txt");
		convertTruth2UTM("results/" + outputPrefix + ".txt", "results/" + outputPrefix + ".txt", dir, KML_FILE);
	}
	
	if (scoreMode) {
		system(("java -Xmx1500m -jar visualizerx.jar -kml training/Challenge1.kml -truth training/ground-truth/Challenge1_Lidar.xyz -novis -solution " + input).c_str());
	}
	
	if (randomMatchMode) {
		VS files = fileList(NITF_DIRECTORY);
		sortNITFFiles(files);
		for (string &s : files) s = s.substr(0, s.S - 4);
		int ID = 0;
		int diff = 1;
		int a = -1;
		while (true) {
			a++;
			int b = a + diff;
			if (b >= files.S) {
				a = -1;
				diff++;
				continue;
			}
			
			if (readRPCOffset(0, files[a]).X == -1) continue;
			if (readRPCOffset(0, files[b]).X == -1) continue;
			system(("./x -r 0 tmp " + files[a] + " " + files[b] + " +" + options + " =" + parameters).c_str());
			system(("./x -v results/tmp.txt images/tmp" + i2s(ID) + ".bmp").c_str());
			system(("echo ID: " + i2s(ID) + " >> " + output).c_str());
			system(("echo " + files[a] + " >> " + output).c_str());
			system(("echo " + files[b] + " >> " + output).c_str());
			system(("./x -s results/tmp.txt >> " + output).c_str());
			system(("echo >> " + output).c_str());
			ID++;
		}
	}
	return 0;
}

