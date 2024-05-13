using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace РГР_ЧМІ__6_Мілса_Сімпсона
{
    internal class Program
    {
        static void Main(string[] args)
        {
            ScottPlot.Plot myPlot = new ScottPlot.Plot();

            Func<double, double, double> f = (t, y) => y + 3 * t - t * t;
            Func<double, double, double> fy = (t, y) => 1;
            ABM abm = new ABM(f, fy, 0, 5);

            abm.WithoutLimitingStep_WithControlParameter(out List<double> X2, out List<double> Y2);
            abm.WithLimitingStep_WithControlParameter(out List<double> X4, out List<double> Y4);
            Func<double, double> fAccuracy = (t) => 2 * Math.Exp(t) + t * t - t - 1;
            Accuracy accuracy = new Accuracy(fAccuracy, abm.a, abm.b, abm.x0, abm.y0, 0.001);
            accuracy.SolveAccuracy(out List<double> X3, out List<double> Y3);

            var line = myPlot.Add.Scatter(X2, Y2);
            var line2 = myPlot.Add.Scatter(X4, Y4);
            var line3 = myPlot.Add.Scatter(X3, Y3);
            line.LineWidth = 4;
            line2.LineWidth = 2;
            line3.LineWidth = 1;
            myPlot.SavePng("mils-simpson.png", 1900, 1800);
        }
        class ABM
        {
            Func<double, double, double> f;
            Func<double, double, double> fy;
            public double x0 = 0;
            public double y0 = 1;

            public double a;
            public double b;
            public int N;
            public double h;

            public ABM(Func<double, double, double> f, Func<double, double, double> fy, double a, double b)
            {
                this.f = f;
                this.fy = fy;
                h = 0.1 / Math.Abs(fy(x0, y0));
                N = (int)((b - a) / h);
                this.a = a;
                this.b = b;
            }
            public void WithoutLimitingStep_WithControlParameter(out List<double> X1, out List<double> Y1)
            {
                ArrayToList(out List<double> X, out List<double> Y);

                int k = 3;
                for (; k < N; k++)
                {
                    X.Add(X[k] + h);
                    double pk = Y[k - 3] + (4 * h / 3) *
                        (
                        2 * f(X[k - 2], Y[k - 2])
                        - f(X[k - 1], Y[k - 1])
                        + 2 * f(X[k], Y[k])
                        );
                    Y.Add
                        (
                            Y[k-1] + (h / 3)
                            * (f(X[k-1], Y[k-1])
                            + 4 * f(X[k], Y[k])
                            + f(X[k+1], pk))
                        );
                }

                X1 = X;
                Y1 = Y;
            }//без управ. парам. з гранич. кроком
            public void WithLimitingStep_WithControlParameter(out List<double> X1, out List<double> Y1)
            {
                ArrayToList(out List<double> X, out List<double> Y);

                int k = 3;
                for (; k < N; k++)
                {
                    X.Add(X[k] + h);
                    double pk1 = Y[k - 3] + (4 * h / 3) *
                        (
                        2 * f(X[k - 2], Y[k - 2])
                        - f(X[k - 1], Y[k - 1])
                        + 2 * f(X[k], Y[k])
                        );
                    double m;
                    if (k == 4)
                        m = pk1 + 9 * (Y[k] - Y[k]) / 121;
                    else
                        m = pk1 + 9 * (Y[k] - pk1) / 121;

                    Y.Add
                        (
                            Y[k - 1] + (h / 3)
                            * (f(X[k - 1], Y[k - 1])
                            + 4 * f(X[k], Y[k])
                            + f(X[k + 1], pk1))
                        );
                }

                X1 = X;
                Y1 = Y;
            }


            private void ArrayToList(out List<double> X, out List<double> Y)
            {
                Hemming hemming = new Hemming(f, h, a, b, x0, y0);
                hemming.Solve(out double[] X1, out double[] Y1);
                X = new List<double>();
                Y = new List<double>();
                for (int i = 0; i < 4; i++)
                {
                    X.Add(X1[i]);
                    Y.Add(Y1[i]);
                }
            }
            class Hemming
            {
                static double eps = 0.00001;

                static Func<double, double, double> f1;
                static Func<double, double, double> f2 = (t, y) => f1(t + eps, y) + f1(t, y + eps) * f1(t, y);
                static Func<double, double, double> f3 = (t, y) => (f1(t + 2 * eps, y) - 2 * f1(t + eps, y)
                + 2 * f1(t - eps, y) - f1(t - 2 * eps, y)) / (2 * eps * eps * eps);

                double x0;
                double y0;

                static double h;

                double[] Xi;
                double[] Yi;
                public Hemming(Func<double, double, double> f1, double h, double a, double b, double x0, double y0)
                {
                    Hemming.h = h * 2;
                    Hemming.f1 = f1;
                    Xi = new double[(int)((b - a) / h)];
                    Yi = new double[(int)((b - a) / h)];
                    this.x0 = x0;
                    this.y0 = y0;
                }
                public void Solve(out double[] Xi1, out double[] Yi1)
                {
                    double yk = y0;

                    Xi[0] = x0;
                    Yi[0] = y0;
                    for (int i = 2; i <= Xi.GetLength(0) / 2; i++)
                    {
                        yk = yk + h * f1(x0, yk) + ((h * h * f2(x0, yk)) / 2) + ((h * h * h * f3(x0, yk) / 6));
                        x0 += h;
                        Xi[i - 1] = x0;
                        Yi[i - 1] = yk;
                    }
                    Xi1 = Xi;
                    Yi1 = Yi;
                }
            }
        }
        class Accuracy
        {
            double a, b;
            double h;
            double x0, y0;
            Func<double, double> f;

            public Accuracy(Func<double, double> f, double a, double b, double x0, double y0, double h)
            {
                this.a = a;
                this.b = b;
                this.h = h * 2;
                this.f = f;
                this.x0 = x0;
                this.y0 = y0;
            }

            public void SolveAccuracy(out List<double> X, out List<double> Y)
            {
                double c = a;
                X = new List<double>();
                Y = new List<double>();
                X.Add(x0);
                Y.Add(y0);
                for (int i = 1; c < b; c += h, i++)
                {
                    X.Add(X[i - 1] + h);
                    Y.Add(f(X[i]));
                }
            }
        }
    }
}
