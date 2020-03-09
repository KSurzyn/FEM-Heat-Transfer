using System;
using Simulation.Utils;

namespace Simulation.FEM {

    public class UniversalElement {
        
        public IntegrationPoint[] Points;
        public int PointsCount;
        public static IntegrationPoint IntegrationPointPC;
        // [ indeks punktu , funkcja ksztaltu ]
        public double[,] NValuesMatrix, dNdKsiMatrix, dNdEtaMatrix;


        public UniversalElement(IntegrationPoint[] points) {
            Points = points;
            PointsCount = points.Length;

            NValuesMatrix = new double[PointsCount, 4];
            dNdKsiMatrix = new double[PointsCount, 4];
            dNdEtaMatrix = new double[PointsCount, 4];
            Calculate();
        }

        private void Calculate() {
            for (int i = 0; i < PointsCount; i++) {
                var ksi = Points[i].Ksi;
                var eta = Points[i].Eta;

                NValuesMatrix[i, 0] = 0.25 * (1 - ksi) * (1 - eta);
                NValuesMatrix[i, 1] = 0.25 * (1 + ksi) * (1 - eta);
                NValuesMatrix[i, 2] = 0.25 * (1 + ksi) * (1 + eta);
                NValuesMatrix[i, 3] = 0.25 * (1 - ksi) * (1 + eta);

                dNdKsiMatrix[i, 0] = -0.25 * (1 - eta);
                dNdKsiMatrix[i, 1] = 0.25 * (1 - eta);
                dNdKsiMatrix[i, 2] = 0.25 * (1 + eta);
                dNdKsiMatrix[i, 3] = -0.25 * (1 + eta);

                dNdEtaMatrix[i, 0] = -0.25 * (1 - ksi);
                dNdEtaMatrix[i, 1] = -0.25 * (1 + ksi);
                dNdEtaMatrix[i, 2] = 0.25 * (1 + ksi);
                dNdEtaMatrix[i, 3] = 0.25 * (1 - ksi);
            }
            
        }
        
      
        public static UniversalElement Create2Point() {
            double p1 = -1d / Math.Sqrt(3);
            double p2 = 1d / Math.Sqrt(3);
            var element = new UniversalElement(
                new[] {
                    new IntegrationPoint(p1, p1, 1, 1), //lewy dolny
                    new IntegrationPoint(p2, p1, 1, 1), //prawy dolny
                    new IntegrationPoint(p2, p2, 1, 1), //prawy gorny
                    new IntegrationPoint(p1, p2, 1, 1)  //lewy gorny
                }
            );
            return element;
        }

        public static double[,] NVector(IntegrationPoint integrationPoint)  // {N} z poszczegilnych N1..N4
        {
            double[,] NVect = new double[4, 1];
            NVect[0, 0] = N1_local(integrationPoint);
            NVect[1, 0] = N2_local(integrationPoint);
            NVect[2, 0] = N3_local(integrationPoint);
            NVect[3, 0] = N4_local(integrationPoint);

            return NVect;
        }

        public static double[,] NVectorTransposed(IntegrationPoint integrationPoint)  // {N} transponowanie
        {
            double[,] NVectTransposed = new double[1, 4];
            NVectTransposed[0, 0] = N1_local(integrationPoint);
            NVectTransposed[0, 1] = N2_local(integrationPoint);
            NVectTransposed[0, 2] = N3_local(integrationPoint);
            NVectTransposed[0, 3] = N4_local(integrationPoint);

            return NVectTransposed;
        }

        public static double N1_local(IntegrationPoint integrationPoint)
        {
            double N1 = 0.25 * (1 - integrationPoint.Ksi) * (1 - integrationPoint.Eta);
            return N1;
        }

        public static double N2_local(IntegrationPoint integrationPoint)
        {
            double N2 = 0.25 * (1 + integrationPoint.Ksi) * (1 - integrationPoint.Eta);
            return N2;
        }

        public static double N3_local(IntegrationPoint integrationPoint)
        {
            double N3 = 0.25 * (1 + integrationPoint.Ksi) * (1 + integrationPoint.Eta);
            return N3;
        }

        public static double N4_local(IntegrationPoint integrationPoint)
        {
            double N4 = 0.25 * (1 - integrationPoint.Ksi) * (1 + integrationPoint.Eta);
            return N4;
        }

      
    }
}