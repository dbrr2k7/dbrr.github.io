---
layout: custompost
title: Phương trình Diophantine tuyến tính
mathjax: true
date: 2023-03-25 21:35 +0700
---

Nguồn: [CP ALGORITHMS](https://cp-algorithms.com/algebra/linear-diophantine-equation.html)

Một phương trình Diophantine tuyến tính (còn gọi là phương trình bậc nhất hai ẩn) là phương trình có dạng:
$ax+by=c$

với $a,b,c$ là các hệ số được cho trước, $x,y$ là các giá trị chưa biết.

Trong bài viết này, chúng ta sẽ tìm hiểu cách để giải quyết các vấn đề thường gặp với phương trình này:
- Tìm một nghiệm bất kỳ của phương trình;
- Tìm tất cả các nghiệm thỏa mãn;

Với:
- Tìm nghiệm trên một khoảng cho trước;
- Tìm nghiệm với giá trị nhỏ nhất của $x+y$.

thì mình sẽ để ở một bài viết khác.

## Trường hợp suy biến
Một trường hợp suy biến mà ta cần quan tâm là trường hợp $a=b=0$. Rõ ràng ta thấy trong trường hợp này phương trình này sẽ vô nghiệm (nếu $c\ne 0$) hoặc vô số nghiệm (nếu $c=0$). Trong phần còn lại của bài viết, ta sẽ bỏ qua trường hợp này.

## Giải kiểu toán
Khi $a\ne 0$ và $b\ne 0$, phương trình $ax+by=c$ có thể được viết lại như sau:
\begin{gather}
ax \equiv c \pmod b,\newline
by \equiv c \pmod a.
\end{gather}

Không mất tính tổng quát, giả sử với $b=0$. Xét phương trình đầu tiên. Khi $a$ và $b$ nguyên tố cùng nhau, nghiệm của phương trình được tính như sau:
$x\equiv ca^{-1} \pmod b$
trong đó $a^{-1}$ là nghịch đảo modulo $b$ của $a$

Khi $a$ và $b$ không nguyên tố cùng nhau (tức là $\gcd(a,b)\ne 1$), giá trị của $ax$ modulo $b$ với mọi số nguyên $x$ sẽ chia hết cho $g=\gcd(a,b)$, vì vậy nghiệm chỉ tồn tại khi $c$ chia hết cho $g$. Trong trường hợp này, một trong các nghiệm của phương trình có thể tìm được bằng cách rút gọn phương trình cho $g$:
$\left(\frac{a}{g}\right)x\equiv \left(\frac{c}{g}\right) \left(\mod \frac{b}{g}\right)$

Vì $g=\gcd(a,b)$ nên $\frac{a}{g}$ và $\frac{b}{g}$ nguyên tố cùng nhau, vì vậy nghiệm của phương trình sẽ có dạng như sau:
\begin{cases}
x \equiv \left(\frac{c}{g}\right)\left(\frac{a}{g}\right)^{-1}\left(\mod{\frac{b}{g}}\right),\\
y = \frac{c-ax}{b}.
\end{cases}

## Giải bằng thuật toán
Để tìm **1 nghiệm bất kỳ** của phương trình Diophantine với 2 ẩn, ta có thể dùng [thuật toán Euclid mở rộng](https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm) (bạn có thể đọc bản tiếng Việt tại [đây](https://vi.wikipedia.org/wiki/Gi%E1%BA%A3i_thu%E1%BA%ADt_Euclid_m%E1%BB%9F_r%E1%BB%99ng), nó cũng được viết khá chi tiết về giải thích và cài đặt). Đầu tiên, giả sử rằng $a$ và $b$ đều không âm. Khi ta áp dụng thuật toán Euclid mở rộng vào $a$ và $b$, ta có thể tìm được ƯCLN của nó và hai số nguyên $x_g$ và $y_g$ sao cho:
$ax_g+by_g=g$

Nếu $c$ chia hết cho $g=\gcd(a,b)$, thì phương trình được cho có một nghiệm, ngược lại phương trình vô nghiệm. Chứng minh điều này khá đơn giản: [tổ hợp tuyến tính](https://vi.wikipedia.org/wiki/T%E1%BB%95_h%E1%BB%A3p_tuy%E1%BA%BFn_t%C3%ADnh) của hai số luôn chia hết cho ƯCLN của chúng.

Bây giờ, giả sử rằng $c$ chia hết cho $g$, ta có:
$a\times x_g\times \frac{c}{g} + b\times y_g \times \frac{c}{g} = c$

Vì thế một trong các nghiệm tồn tại của phương trình sẽ có dạng:
\begin{gather}
x_0 = x_g \times \frac{c}{g}, \newline
y_0 = y_g \times \frac{c}{g}.
\end{gather}

Ý tưởng trên vẫn đúng trong trường hợp $a$ và (hoặc) $b$ là số âm. Ta chỉ cần đổi dấu của chúng trong trường hợp cần thiết.

Ta có thể cài đặt code dựa theo ý tưởng trên như sau: (lưu ý code này không xét trường hợp $a$ và $b$ cùng bằng $0$)

```cpp
int gcd(int a, int b, int& x, int& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    int x1, y1;
    int d = gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return d;
}

bool find_any_solution(int a, int b, int c, int &x0, int &y0, int &g) {
    g = gcd(abs(a), abs(b), x0, y0);
    if (c % g) {
        return false;
    }

    x0 *= c / g;
    y0 *= c / g;
    if (a < 0) x0 = -x0;
    if (b < 0) y0 = -y0;
    return true;
}
```

## Tìm tất cả các nghiệm của phương trình
Từ một nghiệm $\left(x_0,y_0\right)$, ta có thể suy ra được tất cả các nghiệm của phương trình đã cho.

Đặt $g=\gcd(a,b)$, gọi $x_0,y_0$ là các số nguyên thỏa mãn:
$a\times x_0 + b\times y_0=c$

Bây giờ ta thấy rằng khi cộng thêm $\frac{b}{g}$ vào $x_0$ và trừ $y_0$ đi $\frac{a}{g}$ thì phương trình sẽ không đổi:
$a \times \left(x_0 + \frac{b}{g}\right) + b \times \left(y_0 - \frac{a}{g}\right) = a \times x_0 + b \times y_0 + a \times \frac{b}{g} - b \times \frac{a}{g} = c$

Rõ ràng ta thấy, quá trình này có thể được lặp lại thêm một lần nữa, vì thế tất cả các số có dạng:
$x = x_0 + k \cdot \frac{b}{g};$
$y = y_0 - k \cdot \frac{a}{g}$

là các nghiệm của phương trình diophantine đã cho.

Hơn nữa, đây là tập tất cả các nghiệm có thể của phương trình đã cho.
