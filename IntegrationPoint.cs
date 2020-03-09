namespace Simulation.FEM {
    
     public class IntegrationPoint {
        public double Ksi, Eta;
        public double WeightKsi, WeightEta;

        public IntegrationPoint(double ksi, double eta, double weightKsi, double weightEta) {
            Ksi = ksi;
            Eta = eta;
            WeightKsi = weightKsi;
            WeightEta = weightEta;
        }
    }
}