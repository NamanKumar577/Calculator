// server.c - Termux localhost calculator server (no math.h)
// Build: clang -O2 -Wall -Wextra server.c -o calc_server
// Run:   ./calc_server
//
// Open UI: http://127.0.0.1:8080/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>

#define PI  3.141592653589793238462643383279502884
#define TWO_PI (2.0*PI)
#define E   2.718281828459045235360287471352662498

// ---------- small helpers (no math.h) ----------
static double my_abs(double x){ return x < 0 ? -x : x; }

static int is_close_int(double x) {
  double r = (x >= 0) ? (double)(long long)(x + 0.5) : (double)(long long)(x - 0.5);
  return my_abs(x - r) < 1e-9;
}

static long long round_to_ll(double x) {
  return (x >= 0) ? (long long)(x + 0.5) : (long long)(x - 0.5);
}

// Reduce to [-PI, PI] without fmod()
static double reduce_angle(double x){
  long long k = (long long)(x / TWO_PI);
  x -= (double)k * TWO_PI;
  while (x > PI) x -= TWO_PI;
  while (x < -PI) x += TWO_PI;
  return x;
}

// ---------- RK4 integrators ----------
typedef double (*gfun_t)(double t);

// RK4 integrate y' = g(t), y(a)=0, return y(b)
static double rk4_integrate_g(gfun_t g, double a, double b, double h){
  if (h <= 0) h = 0.01;
  double y = 0.0;
  double t = a;
  int forward = (b >= a);
  if (!forward) { double tmp=a; a=b; b=tmp; forward=1; } // ensure a<=b, later flip sign

  t = a;
  while (t < b) {
    double step = h;
    if (t + step > b) step = b - t;

    double k1 = step * g(t);
    double k2 = step * g(t + step*0.5);
    double k3 = step * g(t + step*0.5);
    double k4 = step * g(t + step);

    y += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
    t += step;
  }
  // if original was b<a, integral from a to b is negative of from b to a
  if (b < a) y = -y;
  return y;
}

// RK4 integrate scalar ODE y' = f(t,y), from t=a to t=b with initial y0
typedef double (*ffun_t)(double t, double y);

static double rk4_integrate_f(ffun_t f, double a, double b, double y0, double h){
  if (h <= 0) h = 0.01;
  double t = a;
  double y = y0;

  int forward = (b >= a);
  if (!forward) {
    // integrate backwards by stepping negative
    while (t > b) {
      double step = h;
      if (t - step < b) step = t - b;
      step = -step;

      double k1 = step * f(t, y);
      double k2 = step * f(t + step*0.5, y + k1*0.5);
      double k3 = step * f(t + step*0.5, y + k2*0.5);
      double k4 = step * f(t + step,     y + k3);

      y += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
      t += step;
    }
    return y;
  }

  while (t < b) {
    double step = h;
    if (t + step > b) step = b - t;

    double k1 = step * f(t, y);
    double k2 = step * f(t + step*0.5, y + k1*0.5);
    double k3 = step * f(t + step*0.5, y + k2*0.5);
    double k4 = step * f(t + step,     y + k3);

    y += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
    t += step;
  }
  return y;
}

// ---------- scientific functions (no math.h) ----------

// sin/cos via RK4 on system: y' = z, z' = -y
static void sincos_rk4(double x, double *s, double *c){
  // Use symmetry to keep integration short and accurate
  // sin(-x)=-sin(x), cos(-x)=cos(x)
  int neg = 0;
  x = reduce_angle(x);
  if (x < 0) { neg = 1; x = -x; }

  // further map to [0, PI] using sin(PI - x)=sin(x), cos(PI - x)=-cos(x)
  int flip_cos = 0;
  if (x > PI) x = TWO_PI - x; // not needed after reduce_angle, but safe
  if (x > PI/2) { x = PI - x; flip_cos = 1; }

  double h = 0.005; // good default
  double t = 0.0;
  double y = 0.0; // sin
  double z = 1.0; // cos

  while (t < x) {
    double step = h;
    if (t + step > x) step = x - t;

    // k for y, l for z
    double k1 = step * z;
    double l1 = step * (-y);

    double k2 = step * (z + l1*0.5);
    double l2 = step * (-(y + k1*0.5));

    double k3 = step * (z + l2*0.5);
    double l3 = step * (-(y + k2*0.5));

    double k4 = step * (z + l3);
    double l4 = step * (-(y + k3));

    y += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
    z += (l1 + 2*l2 + 2*l3 + l4) / 6.0;

    t += step;
  }

  if (flip_cos) z = -z;
  if (neg) y = -y;

  *s = y;
  *c = z;
}

static double my_sin(double x){ double s,c; sincos_rk4(x,&s,&c); return s; }
static double my_cos(double x){ double s,c; sincos_rk4(x,&s,&c); return c; }
static double my_tan(double x){
  double s,c; sincos_rk4(x,&s,&c);
  if (my_abs(c) < 1e-12) return (s >= 0 ? 1e308 : -1e308);
  return s/c;
}

// exp(x): y' = y, y(0)=1
static double f_exp(double t, double y){ (void)t; return y; }
static double my_exp(double x){
  double h = 0.01;
  return rk4_integrate_f(f_exp, 0.0, x, 1.0, h);
}

// ln(x): y' = 1/t, y(1)=0  => ln(x) = ∫_1^x 1/t dt
static double g_inv_t(double t){ return 1.0 / t; }
static double my_ln(double x){
  if (x <= 0.0) return 0.0/0.0; // NaN
  double h = 0.001;
  if (x == 1.0) return 0.0;
  if (x > 1.0)  return rk4_integrate_g(g_inv_t, 1.0, x, h);
  // for 0<x<1, ln(x) = - ∫_x^1 1/t dt
  return -rk4_integrate_g(g_inv_t, x, 1.0, h);
}

// sqrt(x) via ODE: y' = 1/(2y), y(1)=1  => y(t)=sqrt(t)
static double f_sqrt(double t, double y){ (void)t; return 1.0/(2.0*y); }
static double my_sqrt(double x){
  if (x < 0.0) return 0.0/0.0; // NaN
  if (x == 0.0) return 0.0;
  double h = 0.001;
  if (x == 1.0) return 1.0;
  return rk4_integrate_f(f_sqrt, 1.0, x, 1.0, h);
}

// arctan(x) = ∫_0^x 1/(1+t^2) dt
static double g_atan(double t){ return 1.0 / (1.0 + t*t); }
static double my_atan(double x){
  double h = 0.001;
  if (x >= 0) return rk4_integrate_g(g_atan, 0.0, x, h);
  return -rk4_integrate_g(g_atan, 0.0, -x, h);
}

// arcsin(x) = ∫_0^x 1/sqrt(1-t^2) dt, |x|<=1
static double g_asin(double t){
  double inside = 1.0 - t*t;
  // inside should be >=0 for |t|<=1
  if (inside <= 0.0) return 1e308; // huge slope near 1
  return 1.0 / my_sqrt(inside);
}
static double my_asin(double x){
  if (x < -1.0 || x > 1.0) return 0.0/0.0;
  double h = 0.0005;
  if (x >= 0) return rk4_integrate_g(g_asin, 0.0, x, h);
  return -rk4_integrate_g(g_asin, 0.0, -x, h);
}
static double my_acos(double x){
  return PI/2.0 - my_asin(x);
}

// factorial for nonnegative integers
static double my_fact(double x){
  if (!is_close_int(x)) return 0.0/0.0;
  long long n = round_to_ll(x);
  if (n < 0) return 0.0/0.0;
  double r = 1.0;
  for (long long i=2;i<=n;i++) r *= (double)i;
  return r;
}

// power: a^b via exp(b*ln(a)) for a>0, else allow integer b if a<0
static double my_pow(double a, double b){
  if (a > 0.0) return my_exp(b * my_ln(a));
  if (a == 0.0){
    if (b > 0.0) return 0.0;
    return 0.0/0.0;
  }
  // a<0
  if (!is_close_int(b)) return 0.0/0.0;
  long long n = round_to_ll(b);
  // integer power
  double base = a;
  long long e = n;
  if (e < 0) { base = 1.0/base; e = -e; }
  double res = 1.0;
  while (e){
    if (e & 1) res *= base;
    base *= base;
    e >>= 1;
  }
  return res;
}

// ---------- Expression parsing (Shunting-yard -> RPN) ----------
typedef enum {
  TK_NUM, TK_OP, TK_LP, TK_RP, TK_FUNC, TK_CONST, TK_COMMA, TK_END
} TokType;

typedef struct {
  TokType type;
  double  num;
  char    op[4];     // "+", "-", "*", "/", "^", "u-" , "!"
  char    name[16];  // sin, ln, pi, etc.
} Token;

static int is_name_char(int c){ return isalpha(c) || isdigit(c) || c=='_'; }

static int tokenize(const char *s, Token *out, int max){
  int n=0;
  const char *p=s;
  while (*p && n<max-1){
    while (isspace((unsigned char)*p)) p++;
    if (!*p) break;

    if (isdigit((unsigned char)*p) || *p=='.'){
      char *end=NULL;
      double v = strtod(p, &end);
      if (end==p) return -1;
      out[n++] = (Token){ .type=TK_NUM, .num=v };
      p=end;
      continue;
    }

    if (isalpha((unsigned char)*p)){
      char buf[16]={0};
      int i=0;
      while (*p && is_name_char((unsigned char)*p) && i<15){
        buf[i++] = (char)tolower((unsigned char)*p);
        p++;
      }
      Token t;
      memset(&t,0,sizeof(t));
      // constants
      if (strcmp(buf,"pi")==0 || strcmp(buf,"e")==0){
        t.type = TK_CONST;
        strncpy(t.name, buf, sizeof(t.name)-1);
      } else {
        t.type = TK_FUNC;
        strncpy(t.name, buf, sizeof(t.name)-1);
      }
      out[n++] = t;
      continue;
    }

    // single-char tokens
    char c=*p;
    if (c=='('){ out[n++] = (Token){.type=TK_LP}; p++; continue; }
    if (c==')'){ out[n++] = (Token){.type=TK_RP}; p++; continue; }
    if (c==','){ out[n++] = (Token){.type=TK_COMMA}; p++; continue; }

    if (c=='+'||c=='-'||c=='*'||c=='/'||c=='^'||c=='!'){
      Token t; memset(&t,0,sizeof(t));
      t.type=TK_OP;
      t.op[0]=c; t.op[1]='\0';
      out[n++]=t;
      p++; continue;
    }
    return -2; // unknown char
  }
  out[n++] = (Token){.type=TK_END};
  return n;
}

static int prec(const Token *t){
  if (t->type != TK_OP) return -1;
  if (strcmp(t->op,"u-")==0) return 5;
  if (strcmp(t->op,"!")==0) return 5; // postfix factorial
  if (strcmp(t->op,"^")==0) return 4;
  if (strcmp(t->op,"*")==0 || strcmp(t->op,"/")==0) return 3;
  if (strcmp(t->op,"+")==0 || strcmp(t->op,"-")==0) return 2;
  return 0;
}
static int right_assoc(const Token *t){
  return (t->type==TK_OP && strcmp(t->op,"^")==0) || (t->type==TK_OP && strcmp(t->op,"u-")==0);
}

static int to_rpn(const Token *in, Token *out, int max){
  Token stack[256];
  int sp=0, nout=0;

  TokType prev = TK_END; // to detect unary minus

  for (int i=0; in[i].type != TK_END; i++){
    Token t = in[i];

    if (t.type==TK_NUM || t.type==TK_CONST){
      out[nout++] = t;
      prev = t.type;
      continue;
    }

    if (t.type==TK_FUNC){
      stack[sp++] = t;
      prev = t.type;
      continue;
    }

    if (t.type==TK_COMMA){
      while (sp>0 && stack[sp-1].type != TK_LP){
        out[nout++] = stack[--sp];
      }
      prev = t.type;
      continue;
    }

    if (t.type==TK_OP){
      // unary minus detection
      if (t.op[0]=='-' && (prev==TK_END || prev==TK_LP || prev==TK_OP || prev==TK_COMMA)){
        strcpy(t.op,"u-");
      }
      // factorial is postfix; treat as operator too
      while (sp>0 && stack[sp-1].type==TK_OP){
        Token top = stack[sp-1];
        int p1 = prec(&t), p2 = prec(&top);
        if ((right_assoc(&t) && p1 < p2) || (!right_assoc(&t) && p1 <= p2)){
          out[nout++] = stack[--sp];
        } else break;
      }
      stack[sp++] = t;
      prev = TK_OP;
      continue;
    }

    if (t.type==TK_LP){
      stack[sp++] = t;
      prev = TK_LP;
      continue;
    }

    if (t.type==TK_RP){
      while (sp>0 && stack[sp-1].type != TK_LP){
        out[nout++] = stack[--sp];
      }
      if (sp==0) return -10; // mismatched
      sp--; // pop '('

      // if function on top, pop it too
      if (sp>0 && stack[sp-1].type==TK_FUNC){
        out[nout++] = stack[--sp];
      }
      prev = TK_RP;
      continue;
    }

    return -20;
  }

  while (sp>0){
    if (stack[sp-1].type==TK_LP) return -11;
    out[nout++] = stack[--sp];
  }
  out[nout++] = (Token){.type=TK_END};
  return nout;
}

static int eval_rpn(const Token *rpn, double *result, char *err, size_t errsz){
  double st[256];
  int sp=0;

  for (int i=0; rpn[i].type != TK_END; i++){
    Token t = rpn[i];

    if (t.type==TK_NUM){
      st[sp++] = t.num;
      continue;
    }
    if (t.type==TK_CONST){
      if (strcmp(t.name,"pi")==0) st[sp++] = PI;
      else if (strcmp(t.name,"e")==0) st[sp++] = E;
      else { snprintf(err,errsz,"Unknown constant"); return 0; }
      continue;
    }

    if (t.type==TK_OP){
      if (strcmp(t.op,"u-")==0){
        if (sp<1){ snprintf(err,errsz,"Unary minus needs 1 operand"); return 0; }
        st[sp-1] = -st[sp-1];
        continue;
      }
      if (strcmp(t.op,"!")==0){
        if (sp<1){ snprintf(err,errsz,"! needs 1 operand"); return 0; }
        st[sp-1] = my_fact(st[sp-1]);
        continue;
      }
      if (sp<2){ snprintf(err,errsz,"Operator %s needs 2 operands", t.op); return 0; }
      double b = st[--sp];
      double a = st[--sp];
      double r = 0.0;

      if (strcmp(t.op,"+")==0) r = a + b;
      else if (strcmp(t.op,"-")==0) r = a - b;
      else if (strcmp(t.op,"*")==0) r = a * b;
      else if (strcmp(t.op,"/")==0) r = a / b;
      else if (strcmp(t.op,"^")==0) r = my_pow(a, b);
      else { snprintf(err,errsz,"Unknown operator"); return 0; }

      st[sp++] = r;
      continue;
    }

    if (t.type==TK_FUNC){
      if (sp<1){ snprintf(err,errsz,"Function %s needs 1 arg", t.name); return 0; }
      double x = st[--sp];
      double r;

      if (strcmp(t.name,"sin")==0) r = my_sin(x);
      else if (strcmp(t.name,"cos")==0) r = my_cos(x);
      else if (strcmp(t.name,"tan")==0) r = my_tan(x);
      else if (strcmp(t.name,"asin")==0 || strcmp(t.name,"arcsin")==0) r = my_asin(x);
      else if (strcmp(t.name,"acos")==0 || strcmp(t.name,"arccos")==0) r = my_acos(x);
      else if (strcmp(t.name,"atan")==0 || strcmp(t.name,"arctan")==0) r = my_atan(x);
      else if (strcmp(t.name,"ln")==0) r = my_ln(x);
      else if (strcmp(t.name,"exp")==0) r = my_exp(x);
      else if (strcmp(t.name,"sqrt")==0) r = my_sqrt(x);
      else { snprintf(err,errsz,"Unknown function: %s", t.name); return 0; }

      st[sp++] = r;
      continue;
    }

    snprintf(err,errsz,"Bad token in RPN");
    return 0;
  }

  if (sp != 1){
    snprintf(err,errsz,"Expression error (stack=%d)", sp);
    return 0;
  }
  *result = st[0];
  return 1;
}

static int eval_expr(const char *expr, double *outv, char *err, size_t errsz){
  Token toks[512], rpn[512];
  int nt = tokenize(expr, toks, 512);
  if (nt < 0){ snprintf(err,errsz,"Tokenize error"); return 0; }
  int nr = to_rpn(toks, rpn, 512);
  if (nr < 0){ snprintf(err,errsz,"Parse error"); return 0; }
  return eval_rpn(rpn, outv, err, errsz);
}

// ---------- minimal HTTP + URL decode ----------
static int hexval(int c){
  if ('0'<=c && c<='9') return c-'0';
  if ('a'<=c && c<='f') return 10 + (c-'a');
  if ('A'<=c && c<='F') return 10 + (c-'A');
  return -1;
}

static void url_decode(char *dst, size_t dsz, const char *src){
  size_t j=0;
  for (size_t i=0; src[i] && j+1<dsz; i++){
    if (src[i]=='%'){
      int h1 = hexval((unsigned char)src[i+1]);
      int h2 = hexval((unsigned char)src[i+2]);
      if (h1>=0 && h2>=0){
        dst[j++] = (char)((h1<<4) | h2);
        i += 2;
      } else {
        dst[j++] = src[i];
      }
    } else if (src[i]=='+'){
      dst[j++] = ' ';
    } else {
      dst[j++] = src[i];
    }
  }
  dst[j] = '\0';
}

static void send_all(int fd, const char *s){
  size_t n = strlen(s);
  while (n){
    ssize_t w = write(fd, s, n);
    if (w <= 0) return;
    s += w; n -= (size_t)w;
  }
}

static void respond_text(int fd, const char *ctype, const char *body){
  char hdr[256];
  snprintf(hdr,sizeof(hdr),
    "HTTP/1.1 200 OK\r\n"
    "Content-Type: %s\r\n"
    "Access-Control-Allow-Origin: *\r\n"
    "Connection: close\r\n\r\n", ctype);
  send_all(fd, hdr);
  send_all(fd, body);
}

static void respond_404(int fd){
  send_all(fd,
    "HTTP/1.1 404 Not Found\r\n"
    "Content-Type: text/plain\r\n"
    "Connection: close\r\n\r\n"
    "404 Not Found");
}

// UI files (served inline for simplicity)
static const char *UI_INDEX =
"<!doctype html>\n"
"<html><head><meta name='viewport' content='width=device-width,initial-scale=1' />\n"
"<title>Termux Calculator</title>\n"
"<link rel='stylesheet' href='/style.css'>\n"
"</head><body>\n"
"<div class='wrap'>\n"
"  <div class='screen'>\n"
"    <div id='expr' class='expr'></div>\n"
"    <div id='out' class='out'>0</div>\n"
"  </div>\n"
"  <div class='grid' id='grid'></div>\n"
"</div>\n"
"<script src='/app.js'></script>\n"
"</body></html>\n";

static const char *UI_CSS =
"body{font-family:system-ui;margin:0;background:#0b0f14;color:#e8eef7;display:flex;justify-content:center}\n"
".wrap{max-width:420px;width:100vw;padding:16px;box-sizing:border-box}\n"
".screen{background:#111826;border-radius:16px;padding:14px;margin-bottom:14px;box-shadow:0 6px 20px rgba(0,0,0,.35)}\n"
".expr{opacity:.7;min-height:22px;word-break:break-all}\n"
".out{font-size:34px;font-weight:700;word-break:break-all}\n"
".grid{display:grid;grid-template-columns:repeat(5,1fr);gap:10px}\n"
"button{border:0;border-radius:14px;padding:14px 10px;font-size:16px;font-weight:650;background:#1a2433;color:#e8eef7;box-shadow:0 6px 16px rgba(0,0,0,.25)}\n"
"button.op{background:#23324a}\n"
"button.eq{background:#2f5cff}\n"
"button.danger{background:#b44}\n";

static const char *UI_JS =
"const exprEl=document.getElementById('expr');\n"
"const outEl=document.getElementById('out');\n"
"const grid=document.getElementById('grid');\n"
"let expr='';\n"
"const btns=[\n"
" ['AC','danger'],['DEL','danger'],['(', 'op'],[')','op'],['^','op'],\n"
" ['sin(','op'],['cos(','op'],['tan(','op'],['ln(','op'],['sqrt(','op'],\n"
" ['asin(','op'],['acos(','op'],['atan(','op'],['exp(','op'],['!','op'],\n"
" ['7'],['8'],['9'],['/','op'],['pi','op'],\n"
" ['4'],['5'],['6'],['*','op'],['e','op'],\n"
" ['1'],['2'],['3'],['-','op'],['.','op'],\n"
" ['0'],['+','op'],['=','eq']\n"
"];\n"
"\n"
"function render(){ exprEl.textContent=expr; }\n"
"function setOut(v){ outEl.textContent=v; }\n"
"\n"
"btns.forEach(([label,cls])=>{\n"
"  const b=document.createElement('button');\n"
"  b.textContent=label;\n"
"  if(cls) b.className=cls;\n"
"  b.onclick=()=>onPress(label);\n"
"  grid.appendChild(b);\n"
"});\n"
"\n"
"async function onPress(k){\n"
"  if(k==='AC'){ expr=''; render(); setOut('0'); return; }\n"
"  if(k==='DEL'){ expr=expr.slice(0,-1); render(); return; }\n"
"  if(k==='='){ await compute(); return; }\n"
"  if(k==='pi') k='pi';\n"
"  if(k==='e') k='e';\n"
"  expr += k;\n"
"  render();\n"
"}\n"
"\n"
"async function compute(){\n"
"  const q=encodeURIComponent(expr);\n"
"  try{\n"
"    const r=await fetch(`/eval?expr=${q}`);\n"
"    const j=await r.json();\n"
"    if(j.ok){ setOut(String(j.value)); }\n"
"    else{ setOut('ERR'); }\n"
"  }catch(e){ setOut('ERR'); }\n"
"}\n"
"\n"
"render();\n";

static void handle_client(int fd){
  char buf[4096];
  ssize_t n = read(fd, buf, sizeof(buf)-1);
  if (n <= 0) return;
  buf[n] = '\0';

  // parse first line: GET /path?query HTTP/1.1
  char method[8]={0}, url[2048]={0};
  if (sscanf(buf, "%7s %2047s", method, url) != 2) return;
  if (strcmp(method,"GET") != 0) return;

  // route
  if (strcmp(url,"/")==0){
    respond_text(fd, "text/html; charset=utf-8", UI_INDEX);
    return;
  }
  if (strcmp(url,"/style.css")==0){
    respond_text(fd, "text/css; charset=utf-8", UI_CSS);
    return;
  }
  if (strcmp(url,"/app.js")==0){
    respond_text(fd, "application/javascript; charset=utf-8", UI_JS);
    return;
  }

  if (strncmp(url,"/eval",5)==0){
    const char *q = strchr(url,'?');
    char expr_raw[2048]={0};
    char expr[2048]={0};

    if (q){
      // find expr=
      const char *p = strstr(q, "expr=");
      if (p){
        p += 5;
        strncpy(expr_raw, p, sizeof(expr_raw)-1);
        // trim other params
        char *amp = strchr(expr_raw,'&');
        if (amp) *amp = '\0';
      }
    }
    url_decode(expr, sizeof(expr), expr_raw);

    double val=0.0;
    char err[128]={0};
    int ok = eval_expr(expr, &val, err, sizeof(err));

    char body[512];
    if (ok){
      snprintf(body,sizeof(body), "{\"ok\":1,\"value\":%.15g}\n", val);
    } else {
      snprintf(body,sizeof(body), "{\"ok\":0,\"error\":\"%s\"}\n", err[0]?err:"bad expression");
    }
    respond_text(fd, "application/json; charset=utf-8", body);
    return;
  }

  respond_404(fd);
}

int main(void){
  int s = socket(AF_INET, SOCK_STREAM, 0);
  if (s < 0){ perror("socket"); return 1; }

  int yes=1;
  setsockopt(s, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(yes));

  struct sockaddr_in addr;
  memset(&addr,0,sizeof(addr));
  addr.sin_family = AF_INET;
  addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK); // 127.0.0.1
  addr.sin_port = htons(8080);

  if (bind(s, (struct sockaddr*)&addr, sizeof(addr)) < 0){
    perror("bind"); return 1;
  }
  if (listen(s, 16) < 0){
    perror("listen"); return 1;
  }

  printf("Calculator server running at http://127.0.0.1:8080/\n");
  fflush(stdout);

  while (1){
    int c = accept(s, NULL, NULL);
    if (c < 0) continue;
    handle_client(c);
    close(c);
  }
  return 0;
}
