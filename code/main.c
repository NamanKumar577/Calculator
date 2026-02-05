#ifndef F_CPU
#define F_CPU 16000000UL
#endif

#include <avr/io.h>
#include <util/delay.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_SIZE 100
#define PI 3.14159265358979323846
#define H  0.01

// ---------------- LCD (Arduino D2..D7 -> PD2..PD7) ----------------
#define LCD_RS PD2
#define LCD_EN PD3
#define LCD_D4 PD4
#define LCD_D5 PD5
#define LCD_D6 PD6
#define LCD_D7 PD7

// ---------------- Keypad pins (AS YOU ASKED) ---------------------
// ROWS: D8..D12  => PB0..PB4
#define ROW0_PB0 0   // D8
#define ROW1_PB1 1   // D9
#define ROW2_PB2 2   // D10
#define ROW3_PB3 3   // D11
#define ROW4_PB4 4   // D12

// COLUMNS: A0..A4 => PC0..PC4
#define COL0_PC0 0   // A0
#define COL1_PC1 1   // A1
#define COL2_PC2 2   // A2
#define COL3_PC3 3   // A3
#define COL4_PC4 4   // A4

static uint8_t mode = 0;
static volatile uint8_t unmatched_brackets = 0;

// ===================== PROTOTYPES (avoid implicit-decl errors) =====================
static void lcd_command(uint8_t cmd);
static void lcd_data(uint8_t data);
static void lcd_init(void);
static void lcd_clear(void);
static void lcd_goto(uint8_t addr);
static void lcd_clear_line(uint8_t line);

static double calc_fast_inv_sqrt(double x);
static double calc_sqrt(double x);
static double calc_sin(double x_target);
static double calc_cos(double x);
static double calc_tan(double x);
static double calc_ln(double x);
static double calc_pow(double x, double w);
static double calc_exp(double x);
static double calc_arctan(double x_target);
static double calc_arcsin(double x_target);
static double calc_arccos(double x);

static int is_func(char c);
static void shunting_yard(const char *expr, char *out);
static double eval_rpn(const char *pfx);

static void render_expr_line0(const char *expr);

static char read_key_for_row(uint8_t row);
static void rows_all_high(void);
static void set_row_low(uint8_t row);
static char keypad_getkey(void);

static void expr_push(char *expr, uint8_t *len, char c);
static void expr_push_func(char *expr, uint8_t *len, char fcode);
static void expr_backspace(char *expr, uint8_t *len);

// ============================ LCD ============================
static void lcd_command(uint8_t cmd) {
  PORTD &= ~(1 << LCD_RS);

  PORTD = (PORTD & 0x0F) | (cmd & 0xF0);
  PORTD |= (1 << LCD_EN); _delay_ms(1);
  PORTD &= ~(1 << LCD_EN);

  PORTD = (PORTD & 0x0F) | ((cmd << 4) & 0xF0);
  PORTD |= (1 << LCD_EN); _delay_ms(1);
  PORTD &= ~(1 << LCD_EN);
}

static void lcd_data(uint8_t data) {
  PORTD |= (1 << LCD_RS);

  PORTD = (PORTD & 0x0F) | (data & 0xF0);
  PORTD |= (1 << LCD_EN); _delay_ms(1);
  PORTD &= ~(1 << LCD_EN);

  PORTD = (PORTD & 0x0F) | ((data << 4) & 0xF0);
  PORTD |= (1 << LCD_EN); _delay_ms(1);
  PORTD &= ~(1 << LCD_EN);
}

static void lcd_init(void) {
  _delay_ms(50);
  lcd_command(0x33);
  lcd_command(0x32);
  lcd_command(0x28);
  lcd_command(0x0C);
  lcd_command(0x06);
  lcd_command(0x01);
  _delay_ms(2);
}

static void lcd_clear(void) {
  lcd_command(0x01);
  _delay_ms(2);
}

static void lcd_goto(uint8_t addr) {
  lcd_command(0x80 | addr);
}

static void lcd_clear_line(uint8_t line) {
  lcd_goto(line ? 0x40 : 0x00);
  for (uint8_t i = 0; i < 16; i++) lcd_data(' ');
  lcd_goto(line ? 0x40 : 0x00);
}

// ============================ MATH ============================
static double calc_fast_inv_sqrt(double x) {
  if (x <= 0.0) return 0.0;

  float y = (float)x;
  float x2 = y * 0.5f;

  union { float f; uint32_t i; } u;
  u.f = y;
  u.i = 0x5f3759df - (u.i >> 1);
  y = u.f;

  y = y * (1.5f - (x2 * y * y));
  return (double)y;
}

static double calc_sqrt(double x) {
  if (x <= 0.0) return 0.0;
  return x * calc_fast_inv_sqrt(x);
}

static double calc_sin(double x_target) {
  double h = 0.01;
  double x = 0.0, y = 0.0, z = 1.0;
  double k1, k2, k3, k4, l1, l2, l3, l4;

  double a = x_target;
  int n = (int)(a / (2.0 * PI));
  a -= n * (2.0 * PI);
  if (a < 0) a += 2.0 * PI;
  x_target = a;

  while (x < x_target) {
    if (x + h > x_target) h = x_target - x;

    k1 = h * z;           l1 = -h * y;
    k2 = h * (z + l1/2);  l2 = -h * (y + k1/2);
    k3 = h * (z + l2/2);  l3 = -h * (y + k2/2);
    k4 = h * (z + l3);    l4 = -h * (y + k3);

    y += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
    z += (l1 + 2*l2 + 2*l3 + l4) / 6.0;
    x += h;
  }
  return y;
}

static double calc_cos(double x) { return calc_sin(PI/2.0 - x); }

static double calc_tan(double x) {
  double c = calc_cos(x);
  if (c == 0.0) return 0.0;
  return calc_sin(x) / c;
}

static double calc_ln(double x) {
  if (x <= 0.0) return 0.0;
  if (x == 1.0) return 0.0;
  if (x < 1.0) return -calc_ln(1.0 / x);

  double x0 = 1.0, y = 0.0;

  while (x0 + H <= x) {
    double k1 = H * (1.0 / x0);
    double k2 = H * (1.0 / (x0 + 0.5*H));
    double k3 = H * (1.0 / (x0 + 0.5*H));
    double k4 = H * (1.0 / (x0 + H));
    y += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
    x0 += H;
  }

  double r = x - x0;
  if (r > 0) {
    double k1 = r * (1.0 / x0);
    double k2 = r * (1.0 / (x0 + 0.5*r));
    double k3 = r * (1.0 / (x0 + 0.5*r));
    double k4 = r * (1.0 / (x0 + r));
    y += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
  }
  return y;
}

static double calc_pow(double x, double w) {
  if (w == 0.0) return 1.0;
  if (x == 0.0) return 0.0;

  double x0 = 1.0, y = 1.0;

  if (x < 1.0) {
    while (x0 - H >= x) {
      double h = -H;
      double k1 = h * (w * y / x0);
      double k2 = h * (w * (y + 0.5*k1) / (x0 + 0.5*h));
      double k3 = h * (w * (y + 0.5*k2) / (x0 + 0.5*h));
      double k4 = h * (w * (y + k3) / (x0 + h));
      y += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
      x0 += h;
    }
    double r = x - x0;
    if (r != 0.0) {
      double h = r;
      double k1 = h * (w * y / x0);
      double k2 = h * (w * (y + 0.5*k1) / (x0 + 0.5*h));
      double k3 = h * (w * (y + 0.5*k2) / (x0 + 0.5*h));
      double k4 = h * (w * (y + k3) / (x0 + h));
      y += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
    }
    return y;
  }

  while (x0 + H <= x) {
    double k1 = H * (w * y / x0);
    double k2 = H * (w * (y + 0.5*k1) / (x0 + 0.5*H));
    double k3 = H * (w * (y + 0.5*k2) / (x0 + 0.5*H));
    double k4 = H * (w * (y + k3) / (x0 + H));
    y += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
    x0 += H;
  }

  double r = x - x0;
  if (r > 0.0) {
    double k1 = r * (w * y / x0);
    double k2 = r * (w * (y + 0.5*k1) / (x0 + 0.5*r));
    double k3 = r * (w * (y + 0.5*k2) / (x0 + 0.5*r));
    double k4 = r * (w * (y + k3) / (x0 + r));
    y += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
  }
  return y;
}

static double calc_exp(double x) {
  const double E = 2.718281828459045;
  return calc_pow(E, x);
}

static double calc_arctan(double x_target) {
  int sign = 1;
  if (x_target < 0.0) { sign = -1; x_target = -x_target; }

  if (x_target > 1.0) {
    double v = calc_arctan(1.0 / x_target);
    return sign * (PI/2.0 - v);
  }

  double x = 0.0, y = 0.0, h = 0.01;
  while (x < x_target) {
    if (x + h > x_target) h = x_target - x;
    double k1 = h / (1.0 + x*x);
    double xm = x + h/2.0;
    double k2 = h / (1.0 + xm*xm);
    double k3 = k2;
    double xe = x + h;
    double k4 = h / (1.0 + xe*xe);
    y += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
    x += h;
  }
  return sign * y;
}

static double calc_arcsin(double x_target) {
  if (x_target < -1.0 || x_target > 1.0) return 0.0;
  if (x_target >= 0.9999) return PI/2.0;
  if (x_target <= -0.9999) return -PI/2.0;

  int sign = 1;
  if (x_target < 0.0) { sign = -1; x_target = -x_target; }

  double x = 0.0, y = 0.0, h = 0.01;
  while (x < x_target) {
    if (x + h > x_target) h = x_target - x;
    double k1 = h * calc_fast_inv_sqrt(1.0 - x*x);
    double xm = x + h/2.0;
    double k2 = h * calc_fast_inv_sqrt(1.0 - xm*xm);
    double k3 = k2;
    double xe = x + h;
    double k4 = h * calc_fast_inv_sqrt(1.0 - xe*xe);
    y += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
    x += h;
  }
  return sign * y;
}

static double calc_arccos(double x) { return (PI/2.0) - calc_arcsin(x); }

// ============================ SHUNTING YARD + RPN ============================
typedef struct { char data[MAX_SIZE]; int top; } OpStack;
static void op_push(OpStack *s, char v) { if (s->top < MAX_SIZE-1) s->data[++s->top] = v; }
static char op_pop(OpStack *s) { return (s->top >= 0) ? s->data[s->top--] : '\0'; }
static char op_peek(OpStack *s) { return (s->top >= 0) ? s->data[s->top] : '\0'; }

static int is_func(char c) { return (c=='s'||c=='c'||c=='t'||c=='e'||c=='l'||c=='z'||c=='y'||c=='x'||c=='q'); }
static int prec(char op) {
  switch(op) {
    case '+': case '-': return 1;
    case '*': case '/': return 2;
    case '^': return 3;
    default: return 0;
  }
}

static void shunting_yard(const char *expr, char *out) {
  OpStack ops; ops.top = -1;
  int j = 0;
  unmatched_brackets = 0;

  for (int i = 0; expr[i] != '\0'; i++) {
    char c = expr[i];

    if (isdigit((unsigned char)c) || c=='.') {
      while (isdigit((unsigned char)expr[i]) || expr[i]=='.') out[j++] = expr[i++];
      out[j++] = ' ';
      i--;
    }
else if (c == 'p') {
  // push pi as a number token into postfix
  const char *ps = "3.141592653589793";
  while (*ps) out[j++] = *ps++;
  out[j++] = ' ';
}
    
    else if (is_func(c)) {
      op_push(&ops, c);
    } else if (c == '(') {
      op_push(&ops, c);
      unmatched_brackets++;
    } else if (c == ')') {
      while (op_peek(&ops) != '(' && ops.top != -1) { out[j++] = op_pop(&ops); out[j++] = ' '; }
      if (op_peek(&ops) == '(') { op_pop(&ops); if (unmatched_brackets) unmatched_brackets--; }
      if (is_func(op_peek(&ops))) { out[j++] = op_pop(&ops); out[j++] = ' '; }
    } else if (strchr("+-*/^", c)) {
      while (ops.top != -1 && prec(op_peek(&ops)) >= prec(c)) { out[j++] = op_pop(&ops); out[j++] = ' '; }
      op_push(&ops, c);
    }
  }

  while (ops.top != -1) { out[j++] = op_pop(&ops); out[j++] = ' '; }
  out[j] = '\0';
}

static double eval_rpn(const char *pfx) {
  double st[MAX_SIZE];
  int top = -1;
  char tok[24];

  while (*pfx) {
    if (*pfx == ' ') { pfx++; continue; }

    if (isdigit((unsigned char)*pfx) || *pfx == '.') {
      int k = 0;
      while (*pfx && *pfx != ' ' && k < 23) tok[k++] = *pfx++;
      tok[k] = '\0';
      st[++top] = atof(tok);
      continue;
    }

    if (strchr("+-*/^", *pfx)) {
      double b = st[top--];
      double a = st[top--];
      switch (*pfx) {
        case '+': st[++top] = a + b; break;
        case '-': st[++top] = a - b; break;
        case '*': st[++top] = a * b; break;
        case '/': st[++top] = (b==0.0)?0.0:(a / b); break;
        case '^': st[++top] = calc_pow(a, b); break;
      }
      pfx++;
      continue;
    }

    if (is_func(*pfx)) {
      double a = st[top--];
      switch (*pfx) {
        case 's': st[++top] = calc_sin(a); break;
        case 'c': st[++top] = calc_cos(a); break;
        case 't': st[++top] = calc_tan(a); break;
        case 'e': st[++top] = calc_exp(a); break;
        case 'l': st[++top] = calc_ln(a); break;
        case 'z': st[++top] = calc_arcsin(a); break;
        case 'y': st[++top] = calc_arccos(a); break;
        case 'x': st[++top] = calc_arctan(a); break;
        case 'q': st[++top] = calc_sqrt(a); break;
      }
      pfx++;
      continue;
    }

    pfx++;
  }

  return (top >= 0) ? st[top] : 0.0;
}

// ============================ DISPLAY ============================
static void render_expr_line0(const char *expr) {
  lcd_clear_line(0);
  lcd_goto(0x00);

  uint8_t col = 0;
  for (uint8_t i = 0; expr[i] && col < 16; i++) {
    char c = expr[i];
if (c == 'p') {          // pi constant
  lcd_data('p');
  if (col < 15) { lcd_data('i'); col++; }
  col++;
  continue;
}
    
    if (is_func(c) && expr[i+1] == '(') {
      const char *s = 0;
      if (c=='s') s="sin(";
      else if (c=='c') s="cos(";
      else if (c=='t') s="tan(";
      else if (c=='e') s="exp(";
      else if (c=='l') s="ln(";
      else if (c=='z') s="asin(";
      else if (c=='y') s="acos(";
      else if (c=='x') s="atan(";
      else if (c=='q') s="sqrt(";

      if (s) {
        while (*s && col < 16) { lcd_data(*s++); col++; }
        i++; // skip '('
        continue;
      }
    }
    lcd_data((uint8_t)c);
    col++;
  }
}

// ============================ KEYPAD (ROWS PB0..PB4, COLS PC0..PC4) ============================
static char read_key_for_row(uint8_t row) {
  const char keys[2][5][5] = {
    { // mode 0
      {'0','1','2','3','4'},
      {'5','6','7','8','9'},
      {'+','-','*','/','s'},
      {'c','t','e','l','C'},
      {'b','.','=','m','p'}
    },
    { // mode 1
      {'0','1','2','3','4'},
      {'5','6','7','8','9'},
      {'(',')','^','f','z'},
      {'y','x','d','g','C'},
      {'b','.','=','m','p'}
    }
  };

  if (!(PINC & (1 << COL0_PC0))) return keys[mode][row][0];
  if (!(PINC & (1 << COL1_PC1))) return keys[mode][row][1];
  if (!(PINC & (1 << COL2_PC2))) return keys[mode][row][2];
  if (!(PINC & (1 << COL3_PC3))) return keys[mode][row][3];
  if (!(PINC & (1 << COL4_PC4))) return keys[mode][row][4];
  return 0;
}

static void rows_all_high(void) {
  PORTB |= (1<<ROW0_PB0)|(1<<ROW1_PB1)|(1<<ROW2_PB2)|(1<<ROW3_PB3)|(1<<ROW4_PB4);
}

static void set_row_low(uint8_t row) {
  rows_all_high();
  PORTB &= ~(1 << row); // row 0..4 => PB0..PB4
}

static char keypad_getkey(void) {
  for (uint8_t row = 0; row < 5; row++) {
    set_row_low(row);
    _delay_ms(2);
    char k = read_key_for_row(row);
    rows_all_high();

    if (k) {
      _delay_ms(25);
      set_row_low(row);
      _delay_ms(2);
      char k2 = read_key_for_row(row);
      rows_all_high();
      if (k2 != k) return 0;

      // wait release
      while (1) {
        uint8_t any = 0;
        for (uint8_t r = 0; r < 5; r++) {
          set_row_low(r);
          _delay_ms(1);
          if (read_key_for_row(r)) any = 1;
          rows_all_high();
          if (any) break;
        }
        if (!any) break;
        _delay_ms(10);
      }
      return k;
    }
  }
  return 0;
}

// ============================ EXPR HELPERS ============================
static void expr_push(char *expr, uint8_t *len, char c) {
  if (*len >= MAX_SIZE-1) return;
  expr[(*len)++] = c;
  expr[*len] = '\0';
}

static void expr_push_func(char *expr, uint8_t *len, char fcode) {
  expr_push(expr, len, fcode);
  expr_push(expr, len, '(');
}

static void expr_backspace(char *expr, uint8_t *len) {
  if (*len == 0) return;
  if (*len >= 2 && expr[*len-1] == '(' && is_func(expr[*len-2])) *len -= 2;
  else *len -= 1;
  expr[*len] = '\0';
}

// ============================ MAIN ============================
int main(void) {
  // LCD pins output
  DDRD |= (1<<LCD_RS)|(1<<LCD_EN)|(1<<LCD_D4)|(1<<LCD_D5)|(1<<LCD_D6)|(1<<LCD_D7);

  // Rows (D8..D12 => PB0..PB4) output HIGH
  DDRB |= (1<<ROW0_PB0)|(1<<ROW1_PB1)|(1<<ROW2_PB2)|(1<<ROW3_PB3)|(1<<ROW4_PB4);
  rows_all_high();

  // Columns (A0..A4 => PC0..PC4) input with pullups
  DDRC &= ~((1<<COL0_PC0)|(1<<COL1_PC1)|(1<<COL2_PC2)|(1<<COL3_PC3)|(1<<COL4_PC4));
  PORTC |=  (1<<COL0_PC0)|(1<<COL1_PC1)|(1<<COL2_PC2)|(1<<COL3_PC3)|(1<<COL4_PC4);

  lcd_init();
  lcd_clear();

  char expr[MAX_SIZE]; expr[0] = '\0';
  uint8_t expr_len = 0;

  render_expr_line0(expr);

  while (1) {
    char key = keypad_getkey();
    if (!key) continue;

    if (key == 'C') {
      expr_len = 0; expr[0] = '\0';
      unmatched_brackets = 0;
      lcd_clear();
      render_expr_line0(expr);
      continue;
    }

    if (key == 'b') {
      expr_backspace(expr, &expr_len);
      render_expr_line0(expr);
      continue;
    }

    if (key == 'm') {
      mode ^= 1;
      lcd_clear_line(1);
      lcd_goto(0x40);
      if (mode) { lcd_data('M'); lcd_data('2'); }
      else      { lcd_data('M'); lcd_data('1'); }
      _delay_ms(250);
      lcd_clear_line(1);
      render_expr_line0(expr);
      continue;
    }

    if (key == '=') {
      if (expr_len == 0) continue;

      char postfix[MAX_SIZE*2];
      memset(postfix, 0, sizeof(postfix));

      shunting_yard(expr, postfix);

      // auto close missing ')'
      if (unmatched_brackets) {
        while (unmatched_brackets && expr_len < MAX_SIZE-1) {
          expr_push(expr, &expr_len, ')');
          unmatched_brackets--;
        }
        render_expr_line0(expr);

        memset(postfix, 0, sizeof(postfix));
        shunting_yard(expr, postfix);
      }

      double ans = eval_rpn(postfix);

      char out[17];
      memset(out, 0, sizeof(out));
      dtostrf(ans, 16, 5, out);

      lcd_clear_line(1);
      lcd_goto(0x40);
      for (uint8_t i = 0; i < 16 && out[i]; i++) lcd_data(out[i]);
      continue;
    }

    // function keys -> add '(' automatically
    if (is_func(key)) {
      expr_push_func(expr, &expr_len, key);
      render_expr_line0(expr);
      continue;
    }

    // normal keys
    expr_push(expr, &expr_len, key);
    render_expr_line0(expr);
  }
}

