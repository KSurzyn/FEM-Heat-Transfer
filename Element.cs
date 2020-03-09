using System;
using Simulation.Utils;

namespace Simulation.FEM {
    public class Element {
        
        public const int ELEMENT_NODES_COUNT = 4;
        
        public int Id;
        public int X, Y;

        public double[,] LocalP;
        public double[,] LocalC;

        public Node[] Nodes;
        public UniversalElement UniversalElement;
 


        public double[][,] Jacobians;
        public double[][,] TransposedJacobians;
        public double[] JacobianDeterminants2D;

        public double[][,] HLocalMatrices;
        public double[,] HLocalMatrix;

        // wektor {N} potrzebny do liczenia macierzy C
        public double[][,] NVector;
        public double[][,] NVectorTransformed;
        public double[][,] MultipliedNVector;  //tablica pomnozonych macierzy
        public IntegrationPoint[][] pc; // w pierwszym wymiarze brzeg a w drugim wymiarze ksi w 0. elemencie i eta w 1. elemencie tablicy

        public Element(int id, int x, int y, Node[] nodes, UniversalElement universalElement) {
            Id = id;
            X = x;
            Y = y;
            Nodes = nodes;
            UniversalElement = universalElement;
            int pcCount = UniversalElement.PointsCount;

            Jacobians = new double[pcCount][,];
            TransposedJacobians = new double[pcCount][,];
            JacobianDeterminants2D = new double[pcCount];

            HLocalMatrices = new double[pcCount][,];
            HLocalMatrix = new double[4, 4];

            for (int i = 0; i < pcCount; i++)
            {
                Jacobians[i] = new double[2, 2];
                TransposedJacobians[i] = new double[2, 2];
                HLocalMatrices[i] = new double[4, 4];
            }
           
        }


        public void Calculate() {
            CalculateJacobians2D();
            TransposedJacobians2D();
            CalculateHLocal(FEM.SimulationData.Instance.Conductivity);
            CalculateC(FEM.SimulationData.Instance.SpecificHeat, FEM.SimulationData.Instance.Ro);
            CalculateH_bc(FEM.SimulationData.Instance.Alfa);
            CalculateP(FEM.SimulationData.Instance.Alfa, FEM.SimulationData.Instance.AmbientTemperature);

        }

        
        private void CalculateJacobians2D() {
            
            for (var i = 0; i < UniversalElement.PointsCount; i++) {
                var jacobian = Jacobians[i];
                //dx/dksi
                jacobian[0, 0] = UniversalElement.dNdKsiMatrix[i, 0] * Nodes[0].GlobalX +
                                 UniversalElement.dNdKsiMatrix[i, 1] * Nodes[1].GlobalX +
                                 UniversalElement.dNdKsiMatrix[i, 2] * Nodes[2].GlobalX +
                                 UniversalElement.dNdKsiMatrix[i, 3] * Nodes[3].GlobalX;
                //dx/deta
                jacobian[1, 0] = UniversalElement.dNdEtaMatrix[i, 0] * Nodes[0].GlobalX +
                                 UniversalElement.dNdEtaMatrix[i, 1] * Nodes[1].GlobalX +
                                 UniversalElement.dNdEtaMatrix[i, 2] * Nodes[2].GlobalX +
                                 UniversalElement.dNdEtaMatrix[i, 3] * Nodes[3].GlobalX;

                //dy/dksi
                jacobian[0, 1] = UniversalElement.dNdKsiMatrix[i, 0] * Nodes[0].GlobalY +
                                 UniversalElement.dNdKsiMatrix[i, 1] * Nodes[1].GlobalY +
                                 UniversalElement.dNdKsiMatrix[i, 2] * Nodes[2].GlobalY +
                                 UniversalElement.dNdKsiMatrix[i, 3] * Nodes[3].GlobalY; 
                
                //dy/deta
                jacobian[1, 1] = UniversalElement.dNdEtaMatrix[i, 0] * Nodes[0].GlobalY +
                                 UniversalElement.dNdEtaMatrix[i, 1] * Nodes[1].GlobalY +
                                 UniversalElement.dNdEtaMatrix[i, 2] * Nodes[2].GlobalY +
                                 UniversalElement.dNdEtaMatrix[i, 3] * Nodes[3].GlobalY;

                JacobianDeterminants2D[i] = jacobian[0, 0] * jacobian[1, 1] - jacobian[0, 1] * jacobian[1, 0];  // wyznacznik jakobianu 2D
         
            }
        }

        private void TransposedJacobians2D() {
            for (int i = 0; i < 4; i++) {
                var jacobian = Jacobians[i];
                TransposedJacobians[i][0, 0] = jacobian[1, 1];
                TransposedJacobians[i][1, 1] = jacobian[0, 0];
                TransposedJacobians[i][0, 1] = -jacobian[0, 1];
                TransposedJacobians[i][1, 0] = -jacobian[1, 0];
            }
        }

        private double DnDx(int pointIndex, int nIndex) {
            //nIndex - numer funkcji ksztaltu
            return 1 / JacobianDeterminants2D[pointIndex] * (TransposedJacobians[pointIndex][0, 0] 
                                                          * UniversalElement.dNdKsiMatrix[pointIndex, nIndex]
                                                          + TransposedJacobians[pointIndex][0, 1] *
                                                          UniversalElement.dNdEtaMatrix[pointIndex, nIndex]);
        }
        
        private double DnDy(int pointIndex, int nIndex) {
            //nIndex - numer funkcji ksztaltu
            return 1 / JacobianDeterminants2D[pointIndex] * (TransposedJacobians[pointIndex][1, 0] 
                                                          * UniversalElement.dNdKsiMatrix[pointIndex, nIndex]
                                                          + TransposedJacobians[pointIndex][1, 1] *
                                                          UniversalElement.dNdEtaMatrix[pointIndex, nIndex]);
        }

        private void CalculateHLocal(double conductivity) {
            for (int i = 0; i < UniversalElement.PointsCount; i++) {
                var dNdX = new[] {
                    DnDx(i, 0),
                    DnDx(i, 1),
                    DnDx(i, 2),
                    DnDx(i, 3)
                };
                var dNdY = new[] {
                    DnDy(i, 0),
                    DnDy(i, 1),
                    DnDy(i, 2),
                    DnDy(i, 3)
                };
                var dNdXMatrix = dNdX.ToMatrix();
                var dNdXMatrixT = dNdX.Transpose();
            
                var dNdYMatrix = dNdY.ToMatrix();
                var dNdYMatrixT = dNdY.Transpose();

                var x = MatrixUtils.MultiplyMatrices(dNdXMatrix, dNdXMatrixT);
                var y = MatrixUtils.MultiplyMatrices(dNdYMatrix, dNdYMatrixT);

                HLocalMatrices[i] = MatrixUtils.AddMatrices(x, y);
            }

            var sum = new double[4, 4];
            for (int i = 0; i < UniversalElement.PointsCount; i++)
            {
                var pc = UniversalElement.Points[i];
                var factor = pc.WeightKsi * pc.WeightEta * JacobianDeterminants2D[i] * conductivity;  
                var matrix = MatrixUtils.MultiplyMatrix(HLocalMatrices[i], factor);
                sum = MatrixUtils.AddMatrices(sum, matrix);  
            }

            HLocalMatrix = sum;
        }

        // c*ro*{N}{N}T  dV
        private void CalculateC(double specificHeat, double ro)  
        {
           NVector = new double[4][,];
           NVectorTransformed = new double[4][,];
           MultipliedNVector = new double[4][,];

            //{N }
            for (int i = 0; i < 4; i++)  
            {
                NVector[i] = new double[4,1]; 
                for (int j = 0; j < 4; j++)  
                {
                    NVector[i][j, 0] = UniversalElement.NValuesMatrix[i, j];   }}
            
            for (int i = 0; i < 4; i++)  
            {
                NVectorTransformed[i] = new double[1,4];
                for (int j = 0; j < 4; j++) 
                {
                    NVectorTransformed[i][0, j] = NVector[i][j, 0]; } }
                
            for (int i = 0; i < 4; i++)
            {
                MultipliedNVector[i] = MatrixUtils.MultiplyMatrix(MatrixUtils.MultiplyMatrices(NVector[i], NVectorTransformed[i]), (specificHeat*ro));
                
            }

            var sum = new double[4, 4];
            for (int i = 0; i < UniversalElement.PointsCount; i++)
            {
                var pc = UniversalElement.Points[i];
                var factor = pc.WeightKsi * pc.WeightEta * JacobianDeterminants2D[i];
                var matrix = MatrixUtils.MultiplyMatrix(MultipliedNVector[i], factor);
                sum = MatrixUtils.AddMatrices(sum, matrix);
            }
            LocalC = sum;

        }

        // alfa*{N}{N}T  dS
        private double[,] CalculateH_bc(int alfa)
        {
           
            pc = new IntegrationPoint[4][]; 
            for (int i = 0; i < 4; i++)
            {
                pc[i] = new IntegrationPoint[2];
            }

                // na kazdy bok przypadają dwa punkty całkowania

                pc[0][0] = new IntegrationPoint(-1 / Math.Sqrt(3), -1, 1, 1);
                pc[0][1] = new IntegrationPoint(1 / Math.Sqrt(3), -1, 1, 1);

                pc[1][0] = new IntegrationPoint(1, -1 / Math.Sqrt(3), 1, 1);
                pc[1][1] = new IntegrationPoint(1, 1 / Math.Sqrt(3), 1, 1);

                pc[2][0] = new IntegrationPoint(1 / Math.Sqrt(3), 1, 1, 1);
                pc[2][1] = new IntegrationPoint(-1 / Math.Sqrt(3), 1, 1, 1);

                pc[3][0] = new IntegrationPoint(-1, 1 / Math.Sqrt(3), 1, 1);
                pc[3][1] = new IntegrationPoint(-1, -1 / Math.Sqrt(3), 1, 1);

       
            double[,] H_bc = new double[4, 4];

            //SPRAWDZANIE WARUNKU BRZEGOWEGO 
            for (int i = 0; i < 4; i++)
            {
                if (Nodes[i].IsBoundary && Nodes[(i + 1) % 4].IsBoundary)
                {
                    Node bc1 = Nodes[i];
                    Node bc2 = Nodes[(i+1)%4];
                    double distance = Math.Sqrt(Math.Pow((bc1.GlobalX - bc2.GlobalX),2) + Math.Pow((bc1.GlobalY - bc2.GlobalY),2));
                    double detJ = distance / 2;

                    IntegrationPoint[] pcForBoundary = pc[i];  
                                                              
                    double[,] matrixFor1pc = new double[4, 4]; 
                    double[,] matrixFor2pc = new double[4, 4];

                    // alfa * {N} * {N}transposed dla kazdego pc
                    matrixFor1pc = MatrixUtils.MultiplyMatrices(UniversalElement.NVector(pcForBoundary[0]), UniversalElement.NVectorTransposed(pcForBoundary[0]));
                    matrixFor2pc = MatrixUtils.MultiplyMatrices(UniversalElement.NVector(pcForBoundary[1]), UniversalElement.NVectorTransposed(pcForBoundary[1]));
                    matrixFor1pc = MatrixUtils.MultiplyMatrix(matrixFor1pc, alfa);
                    matrixFor2pc = MatrixUtils.MultiplyMatrix(matrixFor2pc, alfa);

                    //H_bc += matrixFor1pc * pcForBoundary[0].WeightEta * detJ + matrixFor2pc * pcForBoundary[1].WeightEta * detJ;
                    H_bc = MatrixUtils.AddMatrices(H_bc, MatrixUtils.AddMatrices(MatrixUtils.MultiplyMatrix((MatrixUtils.MultiplyMatrix(matrixFor1pc, pcForBoundary[0].WeightEta)),detJ), MatrixUtils.MultiplyMatrix(MatrixUtils.MultiplyMatrix(matrixFor2pc, pcForBoundary[1].WeightEta),detJ)));
                }
                else
                    continue; 

            }
            this.HLocalMatrix = MatrixUtils.AddMatrices(this.HLocalMatrix, H_bc);
            return H_bc;
        }

        // -alfa*{N} dS
        private void CalculateP(int alfa, double ambientTemperature)
        {

            pc = new IntegrationPoint[4][]; 
            for (int i = 0; i < 4; i++)
            {
                pc[i] = new IntegrationPoint[2];
            }

            // na kazdy bok przypadają dwa punkty całkowania
            pc[0][0] = new IntegrationPoint(-1 / Math.Sqrt(3), -1, 1, 1);
            pc[0][1] = new IntegrationPoint(1 / Math.Sqrt(3), -1, 1, 1);

            pc[1][0] = new IntegrationPoint(1, -1 / Math.Sqrt(3), 1, 1);
            pc[1][1] = new IntegrationPoint(1, 1 / Math.Sqrt(3), 1, 1);

            pc[2][0] = new IntegrationPoint(1 / Math.Sqrt(3), 1, 1, 1);
            pc[2][1] = new IntegrationPoint(-1 / Math.Sqrt(3), 1, 1, 1);

            pc[3][0] = new IntegrationPoint(-1, 1 / Math.Sqrt(3), 1, 1);
            pc[3][1] = new IntegrationPoint(-1, -1 / Math.Sqrt(3), 1, 1);

            double[,] P = new double[4, 1];

            //SPRAWDZANIE WARUNKU BRZEGOWEGO - lewy i prawy dolny róg mają warunek brzegowy ustawiony na 1
            for (int i = 0; i < 4; i++)
            {
                if (Nodes[i].IsBoundary && Nodes[(i + 1) % 4].IsBoundary)
                {
                    Node bc1 = Nodes[i];
                    Node bc2 = Nodes[(i + 1) % 4];
                    double distance = Math.Sqrt(Math.Pow((bc1.GlobalX - bc2.GlobalX), 2) + Math.Pow((bc1.GlobalY - bc2.GlobalY), 2));
                    double detJ = distance / 2;  

                    IntegrationPoint[] pcForBoundary = pc[i];  

                    double[,] matrixFor1pc = new double[4, 1]; 
                    double[,] matrixFor2pc = new double[4, 1];

                    //  temp otoczenia * alfa * {N}  dla kazdego pc
                    matrixFor1pc = UniversalElement.NVector(pcForBoundary[0]);
                    matrixFor2pc = UniversalElement.NVector(pcForBoundary[1]);
                    matrixFor1pc = MatrixUtils.MultiplyMatrix(matrixFor1pc, alfa * ambientTemperature);
                    matrixFor2pc = MatrixUtils.MultiplyMatrix(matrixFor2pc, alfa * ambientTemperature);

                    double[,] integralMatrix = MatrixUtils.AddMatrices(MatrixUtils.MultiplyMatrix(matrixFor1pc, detJ), MatrixUtils.MultiplyMatrix(matrixFor2pc, detJ));
                    //P += matrixFor1pc * pcForBoundary[0].WeightEta * detJ + matrixFor2pc * pcForBoundary[1].WeightEta * detJ;
                    P = MatrixUtils.AddMatrices(P, MatrixUtils.MultiplyMatrix(integralMatrix, -1));
                }
                else
                    continue; 

            }

            LocalP = P; 
        }

    }
}