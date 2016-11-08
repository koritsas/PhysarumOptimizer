package edu.koritsas.slimemold.mstree;

import edu.koritsas.slimemold.shortestpath.PhysarumPolycephalumSP;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.util.FastMath;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;

import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

/**
 * Created by ilias on 8/11/2016.
 */
public abstract class PhysarumPolycephalumMST extends PhysarumPolycephalumSP {
    public PhysarumPolycephalumMST(Graph graph, Node sourceNode, List<Node> sinkNodes, double absoluteThreshold, double relativeThreshold, int numberOfIterations) {
        super(graph, sourceNode, null, absoluteThreshold, relativeThreshold, numberOfIterations);
        this.sinkNodesList=sinkNodes;
    }

    @Override
    public void execute() {
        initializeMaps(graph);

        logger.info("Starting iterations...");

        while (converged==false){

            logger.info("Current iteration: "+iteration);

            Random random = new Random();

            sinkNode =sinkNodesList.get(random.nextInt(sinkNodesList.size()));

            List<Node> allButSinkNodes = getAllButSink(graph);

            RealVector constants = createConstantsVector();

            RealMatrix coefficients = createCoefficientsMatrix();



            // DecompositionSolver solver = new QRDecomposition(coefficients).getSolver();

            DecompositionSolver solver = new SingularValueDecomposition(coefficients).getSolver();

            RealVector solution = solver.solve(constants);

            redefinePressures(allButSinkNodes, solution);

            redefineFlows(graph);

            redefineDiameters(graph);

            converged =checkConvergenceCriterion(previousFluxMap,currentFluxMap,iteration,m_absoluteThreshold,m_relativeThreshold);

            previousFluxMap=new HashMap<>(currentFluxMap);

            iteration++;
        }

    }


}
