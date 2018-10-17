# Matrix

## 拷贝构造函数和移动构造函数

### Intrucion

C++11之前，对象的拷贝控制由三个函数决定：

**拷贝构造函数**（Copy Constructor），

**拷贝赋值运算符**（Copy Assignment Operator）

**析构函数**（Destructor）。

C++11之后，新增加了两个函数：

**移动构造函数**（Move Constructor）

**移动赋值运算符**（Move Assignment Operator）

移动构造函数不可以跨类型移动