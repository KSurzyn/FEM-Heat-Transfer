using Simulation.Utils;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System;

namespace Simulation.FEM
{
    public class Grid
    {
        public Element[] Elements;
        public Node[] Nodes;
        public Grid grid;   
        private GlobalData _globalData;
        public UniversalElement UniversalElement;
        public Grid(GlobalData globalData)
        {
            _globalData = globalData;
            Elements = new Element[globalData.ElementsCount];
            Nodes = new Node[globalData.NodesCount];

            UniversalElement = UniversalElement.Create2Point();
        }

        public void ConstructGrid()
        {
            NodesArray();
            ElementsArray();
            foreach(Element element in Elements)
            {
                element.Calculate();
            }
            agregateCmatrix();
            agregateHmatrix();
            simulation(FEM.SimulationData.Instance.Time, FEM.SimulationData.Instance.DTau, FEM.SimulationData.Instance.Initialtemp);
        }

        public Node GetNodeAt(int x, int y)
        {
            //verify indexes
            int id = x * _globalData.NHeight + y;
            return Nodes[id];
        }

        public Element GetElementAt(int x, int y)
        {
            //verify indexes
            int id = x * (_globalData.NHeight - 1) + y;
            return Elements[id];
        }

        protected void NodesArray()
        {
            int index = 0;
            for (int x = 0; x < _globalData.NLength; x++)
            {
                for (int y = 0; y < _globalData.NHeight; y++)
                {
                    var node = new Node
                    {
                        Id = index,
                        X = x,
                        Y = y,
                        GlobalX = CalculateGlobalX(x),
                        GlobalY = CalculateGlobalY(y),
                        BoundaryCondition = BoundaryCondition(x, y)
                    };
                    Nodes[index] = node;
                    index++;
                }
            }
        }



        protected void ElementsArray()
        {
            int index = 0;
            for (int x = 0; x < _globalData.NLength - 1; x++)
            {
                for (int y = 0; y < _globalData.NHeight - 1; y++)
                {
                    var element = new Element(index, x, y, FindElementNodes(x, y), UniversalElement);
                    Elements[index] = element;
                    index++;
                }
            }
        }
                      
        private Node[] FindElementNodes(int elementX, int elementY)
        {
            var bottomLeft = GetNodeAt(elementX, elementY);
            var topLeft = GetNodeAt(elementX, elementY + 1);
            var topRight = GetNodeAt(elementX + 1, elementY + 1);
            var bottomRight = GetNodeAt(elementX + 1, elementY);
            return new[] { bottomLeft, bottomRight, topRight, topLeft };
        }

        private double CalculateGlobalX(int xIndex)
        {
            return (_globalData.GlobalLengthGrid / (_globalData.NLength - 1)) * xIndex;
        }

        private double CalculateGlobalY(int yIndex)
        {
            return (_globalData.GlobalHeightGrid / (_globalData.NHeight - 1)) * yIndex;
        }

        private int BoundaryCondition(int x, int y)
        {
            if (x == 0 || y == 0
                       || x == _globalData.NLength - 1
                       || y == _globalData.NHeight - 1)
                return 1;
            return 0;
        }


        public double[,] agregateHmatrix()
        {
            double[,] globalH = new double[Nodes.Length, Nodes.Length];

            //petla for po wszystkich wêz³ach siatki
            for (int i = 0; i < Elements.Length; i++)
            {
                    //pobieram id wszystkich wez³ów siatki
                    int[] index_tab = new int[4]{
                        this.Elements[i].Nodes[0].Id,
                        this.Elements[i].Nodes[1].Id,
                        this.Elements[i].Nodes[2].Id,
                        this.Elements[i].Nodes[3].Id
                        };

                    //dodaje poszczegolne elementy do macierzy globalnej H wykorzystuj¹c id wêz³ów
                    for (int y = 0; y < 4; y++)
                    {
                        //petla for po kolumnach macierzy H danego elementu
                        for (int x = 0; x < 4; x++)
                        {
                            //dodajemy element do H zgodnie z id wêz³ów
                            int index_row = index_tab[y];
                            int index_column = index_tab[x];
                            globalH[index_row, index_column] += this.Elements[i].HLocalMatrix[y,x];
                        }
                    }

                
            }
            return globalH;
           
        }


        public double[,] agregateCmatrix()
        {

            double[,] globalC = new double[Nodes.Length, Nodes.Length];

            //petla for po wszystkich wêz³ach siatki
            for (int i = 0; i < Elements.Length; i++)
            {
               
                    //pobieram id wszystkich wez³ów siatki
                    int[] index_tab = new int[4]{
                        this.Elements[i].Nodes[0].Id,
                        this.Elements[i].Nodes[1].Id,
                        this.Elements[i].Nodes[2].Id,
                        this.Elements[i].Nodes[3].Id
                        };


                    //dodaje poszczegolne elementy do macierzy globalnej C wykorzystuj¹c id wêz³ów
                    for (int y = 0; y < 4; y++)
                    {
                        //petla for po kolumnach macierzy C danego elementu
                        for (int x = 0; x < 4; x++)
                        {
                            //dodajemy element do C zgodnie z id wêz³ów
                            int index_row = index_tab[y];
                            int index_column = index_tab[x];
                            globalC[index_row, index_column] += this.Elements[i].LocalC[y, x];
                        }
                    }

            
            }
            return globalC;

        }

        

        public double[,] agregatePvector() 
        {

            double[,] globalP = new double[Nodes.Length, 1];

            //petla for po wszystkich wêz³ach siatki
            for (int i = 0; i < Elements.Length; i++)
            {
               
                    //pobieram id wszystkich wez³ów siatki
                    int[] index_tab = new int[4]{
                        this.Elements[i].Nodes[0].Id,
                        this.Elements[i].Nodes[1].Id,
                        this.Elements[i].Nodes[2].Id,
                        this.Elements[i].Nodes[3].Id
                        };


                    //dodaje poszczegolne elementy do wektora globalnego P wykorzystuj¹c id wêz³ów
                    for (int y = 0; y < 4; y++)
                    {                     
                            //dodajemy element do P zgodnie z id wêz³ów
                            int index_row = index_tab[y];
                            int index_column = 0;
                            globalP[index_row, index_column] += this.Elements[i].LocalP[y, 0];
                        
                    }


            }
            return globalP;

        }


        public void simulation(double time, double dTau, double initialtemp) // (czas symulacji , krok czasowy, temperatura pocz¹tkowa siatki -100)
        {
            //macierze globalne dla przep³ywu nieustalonego ciep³a
            double[,] globalH;
            double[,] globalP;
            double[,] t0 = new double[Nodes.Length, 1];
           
            for (int i = 0; i < t0.Length;  i++)
            {
                t0[i,0] = initialtemp;
            }

            for (double j = dTau; j <= time; j += dTau) // pêtla od czasu dTau do koñca czasu symulacji po kroku dTau
            {
                // [H]+[C]/dT
                globalH = MatrixUtils.AddMatrices(this.agregateHmatrix(), (MatrixUtils.MultiplyMatrix(this.agregateCmatrix(), (1 / dTau))));

                double[,] globalPtemp1 = MatrixUtils.MultiplyMatrix(this.agregateCmatrix(), (1 / dTau));
                double[,] globalPtemp2 = MatrixUtils.MultiplyMatrices(globalPtemp1, t0);
                // {P}+(([C]/dT)*T0)
                globalP = MatrixUtils.AddMatrices(MatrixUtils.MultiplyMatrix(this.agregatePvector(), -1.0), globalPtemp2);
                // odwracanie macierzy globalH za pomoca biblioteki zewnetrznej / [H]{t1} = {P} , szukamy {t1}
                var globalHMath = Matrix<double>.Build.DenseOfArray(globalH);
                var globalPMath = Vector<double>.Build.Dense(globalP.Length, i => globalP[i, 0]);
                var temp1 = globalHMath.Solve(globalPMath);

                // przepisywanie wartoœci z temp1 po ka¿dej iteracji do t0
                for (int i = 0; i < t0.Length; i++)
                {
                    t0[i, 0] = temp1[i];

                }
                // wyswietlanie temperatur
                Console.WriteLine("Time[s]          MinTemp[s]          MaxTemp[s]");
                Console.WriteLine(j + "               " + temp1.Minimum() + "    " + temp1.Maximum());

               
            }

             
        } 


    }
}