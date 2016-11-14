package edu.koritsas.slimemold;

import edu.koritsas.slimemold.shortestpath.PhysarumPolycephalumLangrarianCSP;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Created by ilias on 9/11/2016.
 */
public abstract class TestSlime extends PhysarumPolycephalumLangrarianCSP {
    public TestSlime(Graph graph, Node sourceNode, Node sinkNode, double absoluteThreshold, double relativeThreshold, int numberOfIterations, double λmax, double step) {
        super(graph, sourceNode, sinkNode, absoluteThreshold, relativeThreshold, numberOfIterations, λmax, step);
    }
    public abstract double getConstraintSum(Graph graph);

    @Override
    public void execute() {

        initializeMaps(graph);
        logger.info("Starting iterations...");
        while (λ<=λmax) {
            // for (int i = 0; i <numberOfIterations ; i++) {

            List<Edge> edges = new ArrayList<>(graph.getEdges());
            edges.stream().forEach(edge -> L.put(edge, L.get(edge) + λ * getConstraintSum(getGraph())));
            while (!converged) {




                logger.info("Current iteration: " + iteration);
                List<Node> allNodes = getAllNodes(graph);


                List<Node> allButSinkNodes = getAllButSink(graph);


                RealVector constants = createConstantsVector();

                RealMatrix coefficients = createCoefficientsMatrix();


                // DecompositionSolver solver = new QRDecomposition(coefficients).getSolver();

                DecompositionSolver solver = new SingularValueDecomposition(coefficients).getSolver();

                RealVector solution = solver.solve(constants);

                redefinePressures(allButSinkNodes, solution);

                redefineFlows(graph);

                redefineDiameters(graph);

                converged = checkConvergenceCriterion(previousFluxMap, currentFluxMap, iteration, m_absoluteThreshold, m_relativeThreshold);


                previousFluxMap = new HashMap<>(currentFluxMap);
                iteration++;
                if (iteration==numberOfIterations){
                    break;
                }

            }
            if (iteration==numberOfIterations){
                break;
            }
            converged=false;
            boolean conv =pathViolatesConstraints(getGraph());

            if(conv==false){
                break;
            }



            λ=λ+step;
        }
        System.out.println(converged);

    }

}
