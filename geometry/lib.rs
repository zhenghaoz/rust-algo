// 点
#[derive(Copy,Clone,PartialEq,Debug)]
struct Point {
    x: f64,
    y: f64
}

// 向量
type Vector = Point;

#[derive(PartialEq,Debug)]
enum ClockDirection {
    OnlineBack = 2,
    CounterClockwise = 1,
    Clockwise = -1,
    OnlineFront = -2,
    OnSegment = 0
}

impl Point {

    // 向量加
    fn add(&self, other: &Self) -> Point {
        Point {
            x: self.x + other.x,
            y: self.y + other.y
        }
    }

    // 向量减
    fn sub(&self, other: &Self) -> Point {
        Point {
            x: self.x - other.x,
            y: self.y - other.y
        }
    }

    // 向量长度
    fn abs(&self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    // 标量乘
    fn scale(&self, x: f64) -> Point {
        Point { x: self.x * x, y: self.y * x }
    }

    // 内积
    fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y
    }

    // 外积
    fn cross(&self, other: &Self) -> f64 {
        self.x * other.y - self.y * other.x
    }

    // 计算向量之间余弦
    fn cosine(&self, other: &Self) -> f64 {
        self.dot(other) / self.abs() / other.abs()
    }
    
    // 返回和向量在同一方向的单位向量
    fn unit(&self) -> Point {
        let abs = self.abs();
        Point { x: self.x / abs, y: self.y / abs }
    }

    // 点到点的距离
    fn distance_to_point(&self, other: &Self) -> f64 {
        self.sub(other).abs()
    }

    // 点到直线的距离
    fn distance_to_line(&self, line: &Line) -> f64 {
        self.sub(&line.start).cross(&line.vector()).abs() / line.vector().abs()
    }

    // 点到线段的距离
    fn distance_to_segment(&self, s: &Segment) -> f64 {
        if self.sub(&s.start).dot(&s.end.sub(&s.start)) < std::f64::EPSILON {
            self.distance_to_point(&s.start)
        } else if self.sub(&s.end).dot(&s.start.sub(&s.end)) < std::f64::EPSILON {
            self.distance_to_point(&s.end)
        } else {
            self.distance_to_line(s)
        }
    }

    // 计算向量之间相对方向
    fn ccw(&self, p1: &Self, p2: &Self) -> ClockDirection {
        let v1 = p1.sub(self);
        let v2 = p2.sub(self);
        let c = v1.cross(&v2);
        let d = v1.dot(&v2);
        if c > std::f64::EPSILON {
            ClockDirection::CounterClockwise
        } else if c < -std::f64::EPSILON {
            ClockDirection::Clockwise
        } else if d < -std::f64::EPSILON {
            ClockDirection::OnlineBack
        } else if v1.abs() < v2.abs() {
            ClockDirection::OnlineFront
        } else {
            ClockDirection::OnSegment
        }
    }

    // 计算向量角度
    fn atan2(&self) -> f64 {
        self.y.atan2(self.x)
    }
    
    // 极坐标生成向量
    fn polar(a: f64, r: f64) -> Point {
        Point { x: a.cos() * r, y: a.sin() * r }
    }
}

// 线段
#[derive(Copy,Clone)]
struct Segment {
    start: Point,
    end: Point
}

// 直线
type Line = Segment;

impl Segment {

    // 返回线段对应的向量
    fn vector(&self) -> Vector {
        self.end.sub(&self.start)
    }
    
    // 是否垂直
    fn is_orthogonal(&self, other: &Self) -> bool {
        self.vector().dot(&other.vector()).abs() < std::f64::EPSILON
    }    

    // 是否平行
    fn is_parallel(&self, other: &Self) -> bool {
        self.vector().cross(&other.vector()).abs() < std::f64::EPSILON
    }

    // 计算点p在直线上的投影
    fn projection(&self, p: &Point) -> Point {
        if *p == self.start {
            *p
        } else {
            let hypo = p.sub(&self.start);
            let base = self.vector();
            let t = base.cosine(&hypo) * hypo.abs();
            self.start.add(&self.vector().unit().scale(t))
        }
    }

    // 点p在直线上的映像
    fn reflection(&self, p: &Point) -> Point {
        p.add(&self.projection(&p).sub(&p).scale(2.0))
    }

    // 线段到线段的距离
    fn distance(&self, other: &Self) -> f64 {
        if self.intersect(other) {
            0.0
        } else {
            self.start.distance_to_segment(other)
                .min(self.end.distance_to_segment(other))
                .min(other.start.distance_to_segment(self))
                .min(other.end.distance_to_segment(self))
        }
    }

    // 线段是否相交
    fn intersect(&self, other: &Self) -> bool {
        (self.start.ccw(&self.end, &other.start) as i32) * (self.start.ccw(&self.end, &other.end) as i32) <= 0 &&
        (other.start.ccw(&other.end, &self.start) as i32) * (other.start.ccw(&other.end, &self.end) as i32) <= 0
    }

    // 线段交点
    fn cross_point(&self, s: &Segment) -> Point {
        assert!(self.intersect(s));
        let base = self.vector();
        let d1 = base.cross(&s.start.sub(&self.start)).abs() / base.abs();
        let d2 = base.cross(&s.end.sub(&self.start)).abs() / base.abs();
        let t = d1 / (d1 + d2);
        s.start.add(&s.vector().scale(t))
    }

}

// 圆
#[derive(Copy,Clone)]
struct Circle {
    center: Point,
    radius: f64
}

impl Circle {

    // 圆与直线的交点
    fn cross_points_with_line(&self, line: &Segment) -> (Point, Point) {
        let pr = line.projection(&self.center);
        let e = line.vector().unit();
        let pr_center_length = pr.sub(&self.center).abs();
        let base = (self.radius * self.radius - pr_center_length * pr_center_length).sqrt();
        (pr.add(&e.scale(-base)), pr.add(&e.scale(base)))
    }

    // 圆与圆的交点
    fn cross_points_with_circle(&self, circle: &Self) -> (Point, Point) {
        let d = self.center.sub(&circle.center).abs();
        let a = ((self.radius * self.radius + d * d - circle.radius * circle.radius) / (2.0 * self.radius * d)).acos();
        let t = circle.center.sub(&self.center).atan2();
        (self.center.add(&Point::polar(t + a, self.radius)), self.center.add(&Point::polar(t - a, self.radius)))
    }
}

// 多边形
#[derive(Clone)]
struct Polygen {
    points: Vec<Point>
}

#[derive(PartialEq,Debug)]
enum Containment {
    Out = 0,
    On = 1,
    In = 2
}

impl Polygen {
    // 创建多边形
    fn from(points: Vec<Point>) -> Polygen {
        Polygen { points: points }
    }

    // 创建点集的凸包
    fn convex_hull(_points: &Vec<Point>) -> Polygen {
        let mut points = _points.clone();
        points.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap().then(a.y.partial_cmp(&b.y).unwrap()));
        // 计算上边缘
        let mut hull: Vec<Point> = Vec::new();
        for &point in points.iter() {
            while hull.len() > 1 {
                let p0 = hull[hull.len()-1];
                let p1 = hull[hull.len()-2];
                let direction = p0.ccw(&p1, &point);
                if direction == ClockDirection::CounterClockwise || direction == ClockDirection::OnlineBack {
                    break;
                }
                hull.pop();
            }
            hull.push(point);
        }
        // 计算下边缘
        points.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap().reverse().then(a.y.partial_cmp(&b.y).unwrap().reverse()));
        let n_up_bound = hull.len();
        for &point in points.iter() {
            while hull.len() > n_up_bound {
                let p0 = hull[hull.len()-1];
                let p1 = hull[hull.len()-2];
                let direction = p0.ccw(&p1, &point);
                if direction == ClockDirection::CounterClockwise || direction == ClockDirection::OnlineBack {
                    break;
                }
                hull.pop();
            }
            hull.push(point);
        }
        hull.reverse();
        hull.pop();
        Polygen { points: hull }
    }

    // 计算点是否包含在多边形中
    fn contains(&self, point: &Point) -> Containment {
        let mut cnt = 0;
        for i in 0..self.points.len() {
            let mut a = self.points[i].sub(point);
            let mut b = self.points[(i+1)%self.points.len()].sub(point);
            if a.cross(&b).abs() < std::f64::EPSILON && a.dot(&b) < std::f64::EPSILON {
                return Containment::On;
            }
            if a.y > b.y {
                let temp = a;
                a = b;
                b = temp;
            }
            if a.y < std::f64::EPSILON && b.y > std::f64::EPSILON && a.cross(&b) > std::f64::EPSILON {
                cnt += 1;
            }
        }
        if cnt % 2 == 0 {
            Containment::Out
        } else {
            Containment::In
        }
    }
}

#[cfg(test)]
mod rust_point {

    use {Point,ClockDirection};

    #[test]
    fn test_add() {
        let p1 = Point{ x: 1.0, y: 2.0 };
        let p2 = Point{ x: 3.0, y: 4.0 };
        let ans = p2.add(&p1);
        assert_eq!(4.0, ans.x);
        assert_eq!(6.0, ans.y);
    }

    #[test]
    fn test_sub() {
        let p1 = Point{ x: 1.0, y: 2.0 };
        let p2 = Point{ x: 3.0, y: 4.0 };
        let ans = p2.sub(&p1);
        assert_eq!(2.0, ans.x);
        assert_eq!(2.0, ans.y);
    }

    #[test]
    fn test_abs() {
        let p = Point{ x: 4.0, y: 3.0 };
        assert_eq!(5.0, p.abs());
    }

    #[test]
    fn test_scale() {
        let p = Point{ x: 1.0, y: 2.0 };
        let ans = p.scale(2.0);
        assert_eq!(2.0, ans.x);
        assert_eq!(4.0, ans.y);
    }

    #[test]
    fn test_dot() {
        let p1 = Point{ x: 25.0, y: 0.0 };
        let p2 = Point{ x: 9.0, y: 12.0 };
        assert_eq!(25.0*9.0, p1.dot(&p2));
    }

    #[test]
    fn test_cross() {
        let p1 = Point{ x: 25.0, y: 0.0 };
        let p2 = Point{ x: 9.0, y: 12.0 };
        assert_eq!(20.0*15.0, p1.cross(&p2));
    }

    #[test]
    fn test_ccw() {
        let p0 = Point { x: 0.0, y: 0.0 };
        let p1 = Point { x: 2.0, y: 0.0 };
        assert_eq!(ClockDirection::CounterClockwise, p0.ccw(&p1,&Point{x:-1.0,y:1.0}));
        assert_eq!(ClockDirection::Clockwise, p0.ccw(&p1,&Point{x:-1.0,y:-1.0}));
        assert_eq!(ClockDirection::OnlineBack, p0.ccw(&p1,&Point{x:-1.0,y:0.0}));
        assert_eq!(ClockDirection::OnSegment, p0.ccw(&p1,&Point{x:0.0,y:0.0}));
        assert_eq!(ClockDirection::OnlineFront, p0.ccw(&p1,&Point{x:3.0,y:0.0}));
    }

}

#[cfg(test)]
mod test_segment {

    use {Point, Segment};

    #[test]
    fn test_orthogonal() {
        let s1 = Segment { 
            start: Point { x: 0.0, y: 0.0 },
            end: Point { x: 3.0, y: 0.0 }
        };
        let s2 = Segment { 
            start: Point { x: 1.0, y: 1.0 },
            end: Point { x: 1.0, y: 4.0 }
        };
        assert!(!s1.is_parallel(&s2));
        assert!(s1.is_orthogonal(&s2));
    }

    #[test]
    fn test_parallel() {
        let s1 = Segment { 
            start: Point { x: 0.0, y: 0.0 },
            end: Point { x: 3.0, y: 0.0 }
        };
        let s2 = Segment { 
            start: Point { x: 0.0, y: 0.0 },
            end: Point { x: 3.0, y: 0.0 }
        };
        assert!(s1.is_parallel(&s2));
        assert!(!s1.is_orthogonal(&s2));
    }

    #[test]
    fn test_projection() {
        let s = Segment {
            start: Point { x: 0.0, y: 0.0 },
            end: Point { x: 2.0, y: 0.0 }
        };
        assert_eq!(Point{x:-1.0,y:0.0}, s.projection(&Point{x:-1.0,y:1.0}));
        assert_eq!(Point{x:0.0,y:0.0}, s.projection(&Point{x:0.0,y:1.0}));
        assert_eq!(Point{x:1.0,y:0.0}, s.projection(&Point{x:1.0,y:1.0}));
    }

    #[test]
    fn test_reflection() {
        let s = Segment {
            start: Point { x: 0.0, y: 0.0 },
            end: Point { x: 2.0, y: 0.0 }
        };
        assert_eq!(Point{x:-1.0,y:-1.0}, s.reflection(&Point{x:-1.0,y:1.0}));
        assert_eq!(Point{x:0.0,y:-1.0}, s.reflection(&Point{x:0.0,y:1.0}));
        assert_eq!(Point{x:1.0,y:-1.0}, s.reflection(&Point{x:1.0,y:1.0}));
    }

    #[test]
    fn test_distance() {

        // #1
        let s1a = Segment {
            start: Point { x: 0.0, y: 0.0 },
            end: Point { x: 1.0, y: 0.0 }
        };
        let s1b = Segment {
            start: Point { x: 0.0, y: 1.0 },
            end: Point { x: 1.0, y: 1.0 }
        };
        assert_eq!(1.0, s1a.distance(&s1b));

        // #2
        let s2a = Segment {
            start: Point { x: 0.0, y: 0.0 },
            end: Point { x: 1.0, y: 0.0 }
        };
        let s2b = Segment {
            start: Point { x: 2.0, y: 1.0 },
            end: Point { x: 1.0, y: 2.0 }
        };
        assert!((s2a.distance(&s2b) - 1.4142135624).abs() < 0.00000001);

        // #3
        let s3a = Segment {
            start: Point { x: -1.0, y: 0.0 },
            end: Point { x: 1.0, y: 0.0 }
        };
        let s3b = Segment {
            start: Point { x: 0.0, y: 1.0 },
            end: Point { x: 0.0, y: -1.0 }
        };
        assert_eq!(0.0, s3a.distance(&s3b));
    }

    #[test]
    fn test_intersect() {
        
        // #1
        let s1a = Segment {
            start: Point { x: 0.0, y: 0.0 },
            end: Point { x: 3.0, y: 0.0 }
        };
        let s1b = Segment {
            start: Point { x: 1.0, y: 1.0 },
            end: Point { x: 2.0, y: -1.0 }
        };
        assert!(s1a.intersect(&s1b));

        // #2
        let s2a = Segment {
            start: Point { x: 0.0, y: 0.0 },
            end: Point { x: 3.0, y: 0.0 }
        };
        let s2b = Segment {
            start: Point { x: 3.0, y: 1.0 },
            end: Point { x: 3.0, y: -1.0 }
        };
        assert!(s2a.intersect(&s2b));

        // #3
        let s3a = Segment {
            start: Point { x: 0.0, y: 0.0 },
            end: Point { x: 3.0, y: 0.0 }
        };
        let s3b = Segment {
            start: Point { x: 3.0, y: -2.0 },
            end: Point { x: 5.0, y: 0.0 }
        };
        assert!(!s3a.intersect(&s3b));

    }

    #[test]
    fn test_cross_point() {

        // #1
        let s1a = Segment {
            start: Point { x: 0.0, y: 0.0 },
            end: Point { x: 2.0, y: 0.0 }
        };
        let s1b = Segment {
            start: Point { x: 1.0, y: 1.0 },
            end: Point { x: 1.0, y: -1.0 }
        };
        assert_eq!(Point{x:1.0,y:0.0}, s1a.cross_point(&s1b));

        // #2
        let s2a = Segment {
            start: Point { x: 0.0, y: 0.0 },
            end: Point { x: 1.0, y: 1.0 }
        };
        let s2b = Segment {
            start: Point { x: 0.0, y: 1.0 },
            end: Point { x: 1.0, y: 0.0 }
        };
        assert_eq!(Point{x:0.5,y:0.5}, s2a.cross_point(&s2b));

        // #3
        let s3a = Segment {
            start: Point { x: 0.0, y: 0.0 },
            end: Point { x: 1.0, y: 1.0 }
        };
        let s3b = Segment {
            start: Point { x: 1.0, y: 0.0 },
            end: Point { x: 0.0, y: 1.0 }
        };
        assert_eq!(Point{x:0.5,y:0.5}, s3a.cross_point(&s3b));
    }

}

#[cfg(test)]
mod test_circle {

    use {Point, Segment, Circle};

    #[test]
    fn test_cross_points_with_line() {
        let c = Circle{ center: Point{ x: 2.0, y: 1.0 }, radius: 1.0 };
        let s1 = Segment {
            start: Point { x: 0.0, y: 1.0 },
            end: Point { x: 4.0, y: 1.0 }
        };
        let (p1, p2) = c.cross_points_with_line(&s1);
        assert!((p1.x - 1.00000000).abs() < 0.000001);
        assert!((p1.y - 1.00000000).abs() < 0.000001);
        assert!((p2.x - 3.00000000).abs() < 0.000001);
        assert!((p2.y - 1.00000000).abs() < 0.000001);
        let s2 = Segment {
            start: Point { x: 3.0, y: 0.0 },
            end: Point { x: 3.0, y: 3.0 }
        };
        let (p3, p4) = c.cross_points_with_line(&s2);
        assert!((p3.x - 3.00000000).abs() < 0.000001);
        assert!((p3.y - 1.00000000).abs() < 0.000001);
        assert!((p4.x - 3.00000000).abs() < 0.000001);
        assert!((p4.y - 1.00000000).abs() < 0.000001);
    }

    #[test]
    fn cross_points_with_circle() {
        // case 1
        let c1 = Circle{ center: Point{ x: 0.0, y: 0.0 }, radius: 2.0 };
        let c2 = Circle{ center: Point{ x: 2.0, y: 0.0 }, radius: 2.0 };
        let (p1, p2) = c1.cross_points_with_circle(&c2);
        assert!((p1.x - 1.00000000).abs() < 0.000001);
        assert!((p1.y - 1.73205080).abs() < 0.000001);
        assert!((p2.x - 1.00000000).abs() < 0.000001);
        assert!((p2.y + 1.73205080).abs() < 0.000001);
        // case 2
        let c3 = Circle{ center: Point{ x: 0.0, y: 0.0 }, radius: 2.0 };
        let c4 = Circle{ center: Point{ x: 0.0, y: 3.0 }, radius: 1.0 };
        let (p3, p4) = c3.cross_points_with_circle(&c4);
        assert!((p3.x - 0.00000000).abs() < 0.000001);
        assert!((p3.y - 2.00000000).abs() < 0.000001);
        assert!((p4.x - 0.00000000).abs() < 0.000001);
        assert!((p4.y - 2.00000000).abs() < 0.000001);
    }

}

#[cfg(test)]
mod test_polygen {

    use {Point, Polygen, Containment};

    #[test]
    fn test_convex_hull() {
        // case 1
        let poly1 = Polygen::convex_hull(&vec!(
            Point{x:2.0, y:1.0},
            Point{x:0.0, y:0.0},
            Point{x:1.0, y:2.0},
            Point{x:2.0, y:2.0},
            Point{x:4.0, y:2.0},
            Point{x:1.0, y:3.0},
            Point{x:3.0, y:3.0},
        ));
        assert_eq!(Point{x:0.0,y:0.0}, poly1.points[0]);
        assert_eq!(Point{x:2.0,y:1.0}, poly1.points[1]);
        assert_eq!(Point{x:4.0,y:2.0}, poly1.points[2]);
        assert_eq!(Point{x:3.0,y:3.0}, poly1.points[3]);
        assert_eq!(Point{x:1.0,y:3.0}, poly1.points[4]);

        // case 2
        let poly2 = Polygen::convex_hull(&vec!(
            Point{x:0.0, y:0.0},
            Point{x:2.0, y:2.0},
            Point{x:0.0, y:2.0},
            Point{x:0.0, y:1.0},
        ));
        assert_eq!(Point{x:0.0,y:0.0}, poly2.points[0]);
        assert_eq!(Point{x:2.0,y:2.0}, poly2.points[1]);
        assert_eq!(Point{x:0.0,y:2.0}, poly2.points[2]);
        assert_eq!(Point{x:0.0,y:1.0}, poly2.points[3]);
    }

    #[test]
    fn test_contains() {

        let poly = Polygen::from(vec!(
            Point{x:0.0, y:0.0},
            Point{x:3.0, y:1.0},
            Point{x:2.0, y:3.0},
            Point{x:0.0, y:3.0},
        ));
        assert_eq!(Containment::In, poly.contains(&Point{x:2.0,y:1.0}));
        assert_eq!(Containment::On, poly.contains(&Point{x:0.0,y:2.0}));
        assert_eq!(Containment::Out, poly.contains(&Point{x:3.0,y:2.0}));
    }

}
