package examples;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import org.cloudbus.cloudsim.Cloudlet;
import org.cloudbus.cloudsim.CloudletSchedulerTimeShared;
import org.cloudbus.cloudsim.Datacenter;
import org.cloudbus.cloudsim.DatacenterBroker;
import org.cloudbus.cloudsim.DatacenterCharacteristics;
import org.cloudbus.cloudsim.Host;
import org.cloudbus.cloudsim.Log;
import org.cloudbus.cloudsim.Pe;
import org.cloudbus.cloudsim.Storage;
import org.cloudbus.cloudsim.UtilizationModel;
import org.cloudbus.cloudsim.UtilizationModelFull;
import org.cloudbus.cloudsim.Vm;
import org.cloudbus.cloudsim.VmAllocationPolicy;
import org.cloudbus.cloudsim.VmAllocationPolicySimple;
import org.cloudbus.cloudsim.VmSchedulerTimeShared;
import org.cloudbus.cloudsim.core.CloudSim;
import org.cloudbus.cloudsim.core.CloudSimShutdown;
import org.cloudbus.cloudsim.provisioners.BwProvisionerSimple;
import org.cloudbus.cloudsim.provisioners.PeProvisionerSimple;
import org.cloudbus.cloudsim.provisioners.RamProvisionerSimple;

public class pso_example {
	/** The cloudlet list. */
	private static ArrayList<Cloudlet> cloudletList;

	/** The vmlist. */
	private static List<Vm> vmlist;

	/**
	 * Creates main() to run this example
	 */
	
    public static void main(String[] args) {
    	
    	Log.printLine("Starting PSO_Algorithm...");

        try {
        	// First step: Initialize the CloudSim package. It should be called
        	// before creating any entities.
        	int num_user = 1;   // number of cloud users
        	Calendar calendar = Calendar.getInstance();
        	boolean trace_flag = false;  // mean trace events

			//VM description
			int mips = 1000;
			long size = 20000; //image size (MB)
			int ram = 2048; //vm memory (MB)
			long bw = 1024;
			int pesNumber = 1; //number of cpus
			String vmm = "Xen"; //VMM name
			
			//Fifth step: Create 100 Cloudlets (CREATE AFTER INIT PSO!!!)

            // Step 4: Create a swarm and run the simulation
			// Set the PSO parameters
			int numParticles = 20;
			double inertiaWeight = 0.729;
			double cognitiveCoefficient = 1.49445;
			double socialCoefficient = 1.49445;
			int numCloudlets = 100;
			int numVms = 36;
			double[] velocityLowerBounds = new double[numCloudlets * numVms];
			double[] velocityUpperBounds = new double[numCloudlets * numVms];
			Arrays.fill(velocityLowerBounds, -1.0);
			Arrays.fill(velocityUpperBounds, 1.0);

			// Initialize the particles with random positions and velocities
			int[][] particlePositions = new int[numParticles][numCloudlets];
			int[][] particleVelocities = new int[numParticles][numCloudlets * numVms];
			Random random = new Random();
			for (int i = 0; i < numParticles; i++) {
			    for (int j = 0; j < numCloudlets; j++) {
			        particlePositions[i][j] = random.nextInt(numVms);
			    }
			    for (int j = 0; j < numCloudlets * numVms; j++) {
			        particleVelocities[i][j] = random.nextInt(2) - 1;
			    }
			}
				
			// Run the PSO algorithm for a specified number of iterations
			int numIterations = 50;
			int[] globalBestPosition = new int[numCloudlets];
			double globalBestFitness = Double.MAX_VALUE;
			for (int iteration = 0; iteration < numIterations; iteration++) {
			    for (int i = 0; i < numParticles; i++) {

					// Create Entities
					CloudSim.init(num_user, calendar, trace_flag);
					Datacenter datacenter0 = createDatacenter("Datacenter_0");
					DatacenterBroker broker = createBroker();
					int brokerId = broker.getId();
					vmlist = new ArrayList<Vm>();
					int vmid = 0;
					for (int it = 0; it < 36; it++) {
						Vm vm = new Vm(vmid, brokerId, mips, pesNumber, ram, bw, size, vmm, new CloudletSchedulerTimeShared());
						vmlist.add(vm);
						vmid++;
					}
					broker.submitVmList(vmlist);

					// Initial particle position and velocity
			        int[] particlePosition = particlePositions[i];
			        int[] particleVelocity = particleVelocities[i];

					// Create Cloudlets
					cloudletList = new ArrayList<Cloudlet>();
					createCloudletArray(cloudletList, pesNumber, brokerId);
					// Submit cloudlet list to the broker
					broker.submitCloudletList(cloudletList);

					// Bind Cloudlet
					for (int j = 0; j < numCloudlets; j++) {
						int cloudlet_id = cloudletList.get(j).getCloudletId();
						int vm_id = vmlist.get(particlePosition[j]).getId();
						broker.bindCloudletToVm(cloudlet_id, vm_id);
					}
					
					// Start Simulation 
					CloudSim.startSimulation();

					double particleFitness = calculateFitness(particlePosition, datacenter0, cloudletList, numCloudlets, broker);
					//System.out.println("particle fitness " + particleFitness);
			        int[] personalBestPosition = particlePosition;
			        double personalBestFitness = particleFitness;
			        if (personalBestFitness < globalBestFitness) {
			            globalBestPosition = personalBestPosition;
			            globalBestFitness = personalBestFitness;
			        }
			        for (int j = 0; j < numCloudlets; j++) {
			            double r1 = random.nextDouble();
			            double r2 = random.nextDouble();
			            double cognitiveComponent = cognitiveCoefficient * r1 * (personalBestPosition[j] - particlePosition[j]);
			            double socialComponent = socialCoefficient * r2 * (globalBestPosition[j] - particlePosition[j]);
			            particleVelocity[j] = (int) (inertiaWeight * particleVelocity[j] + cognitiveComponent + socialComponent);
			            particleVelocity[j] = Math.max(particleVelocity[j], (int) velocityLowerBounds[j]);
			            particleVelocity[j] = Math.min(particleVelocity[j], (int) velocityUpperBounds[j]);
						System.out.println("Velocity: " + particleVelocity[j]);
			            particlePosition[j] += particleVelocity[j];
			            particlePosition[j] = Math.max(particlePosition[j], 0);
			            particlePosition[j] = Math.min(particlePosition[j], numVms - 1);
			        }

			    }
			}

			CloudSim.stopSimulation();

			// Simulate PSO
			CloudSim.init(num_user, calendar, trace_flag);
			Datacenter datacenter0 = createDatacenter("Datacenter_0");
			DatacenterBroker broker = createBroker();
			int brokerId = broker.getId();
			vmlist = new ArrayList<Vm>();
			int vmid = 0;
			for (int it = 0; it < 36; it++) {
				Vm vm = new Vm(vmid, brokerId, mips, pesNumber, ram, bw, size, vmm, new CloudletSchedulerTimeShared());
				vmlist.add(vm);
				vmid++;
			}
			broker.submitVmList(vmlist);

			// Create Cloudlets
			cloudletList = new ArrayList<Cloudlet>();
			createCloudletArray(cloudletList, pesNumber, brokerId);
			// Submit cloudlet list to the broker
			broker.submitCloudletList(cloudletList);

			// Assign the cloudlets to the VMs according to the best position found by the PSO algorithm
			for (int i = 0; i < numCloudlets; i++) {
				int cloudlet_id = cloudletList.get(i).getCloudletId();
				int vm_id = vmlist.get(globalBestPosition[i]).getId();
				//System.out.println("Best Position: " + globalBestPosition[i]);
				broker.bindCloudletToVm(cloudlet_id, vm_id);
			}
//
			// Step 5: Start the CloudSim simulation
			CloudSim.startSimulation();
			// Print the results
			List<Cloudlet> finishedCloudlets = new ArrayList<>();

			for (Cloudlet cloudlet : cloudletList) {
			    int cloudletStatus = cloudlet.getCloudletStatus();
			    if (cloudletStatus == Cloudlet.SUCCESS) {
			        finishedCloudlets.add(cloudlet);
			    }
			}

			// Step 6: Stop the CloudSim simulation
			CloudSim.stopSimulation();
			double totalExecutionTime = 0.0;
			for (Cloudlet cloudlet : finishedCloudlets) {
			totalExecutionTime += cloudlet.getActualCPUTime();
			}
			double averageExecutionTime = totalExecutionTime / finishedCloudlets.size();
			System.out.println("Average execution time per cloudlet: " + averageExecutionTime);

        	} catch (Exception e) {
            e.printStackTrace();
            System.err.println("Simulation terminated due to an unexpected error");
        	}
    }	
    
	private static void createCloudletArray(ArrayList<Cloudlet> cloudletList, int pesNumber, int brokerId) {

		// Clear previous cloudlets
		if (cloudletList != null) {
			cloudletList.clear();
		}
		

		int cloudLength[] = {127000, 83000, 63000, 63000, 65000, 15000, 91000, 85000, 101000, 40000, 40000, 712500, 59000, 113000, 87000, 53000, 47000, 127000, 45000, 27500, 115000, 525000, 117000, 105000, 65000, 77000, 47000, 89000, 91000, 105000, 117000, 51000, 15000, 337500, 129000, 47000, 109000, 109000, 135000, 89000, 101000, 95000, 125000, 127000, 83000, 900000, 712500, 83000, 101000, 59000, 105000, 59000, 47000, 125000, 105000, 712500, 45000, 109000, 67000, 61000, 67000, 127000, 87000, 107000, 135000, 101000, 95000, 127000, 525000, 51000, 91000, 91000, 93000, 63000, 51000, 91000, 79000, 77000, 91000, 111000, 85000, 71000, 125000, 712500, 49000, 47000, 40000, 87000, 125000, 75000, 337500, 61000, 93000, 49000, 113000, 115000, 525000, 55000, 65000, 79000};
		
		//Cloudlet properties
		int id = 0;
		long fileSize = 300;
		long outputSize = 300;
		UtilizationModel utilizationModel = new UtilizationModelFull();
		for (int i = 0; i < 100; i++) {
			int length = cloudLength[i];
			Cloudlet cloudlet = new Cloudlet(id, length, pesNumber, fileSize, outputSize, utilizationModel, utilizationModel, utilizationModel);
			cloudlet.setUserId(brokerId);
			cloudletList.add(cloudlet);
			id++;
		}

			//Cloudlet cloudlet1 = new Cloudlet(id, length, pesNumber, fileSize, outputSize, utilizationModel, utilizationModel, utilizationModel);
			//cloudlet1.setUserId(brokerId);

			//id++;
			//Cloudlet cloudlet2 = new Cloudlet(id, length, pesNumber, fileSize, outputSize, utilizationModel, utilizationModel, utilizationModel);
			//cloudlet2.setUserId(brokerId);

			//add the cloudlets to the list
			//cloudletList.add(cloudlet1);
			//cloudletList.add(cloudlet2);

	}

	// Fitness function
    private static double calculateFitness(int[] particlePosition, Datacenter datacenter0, List<Cloudlet> cloudletList, int numCloudlets, DatacenterBroker broker) {
		double result;

    	// Calculate MakeSpan
		double makespan = 0.0;
		for (Cloudlet cloudlet : cloudletList) {
			if (cloudlet.getFinishTime() > makespan) {
				makespan = cloudlet.getFinishTime();
			}
		}
		System.out.println("MakeSpan: " + makespan);

		// Calculate Resource Usage
		List<Host> hostList = datacenter0.getHostList();
		double totalResourceUsage = 0.0;
		//for (Host host : hostList) {
		//	double utilization = host.getTotalUtilizationOfCpu(CloudSim.clock());
		//	totalUtilization += utilization;
		//	numHosts++;
		//}
		double totalUtilization = 0.0;
		for (Host host : hostList) {
			for (Pe pe : host.getPeList()) {
				totalUtilization += pe.getPeProvisioner().getUtilization();
			}
			totalResourceUsage += totalUtilization / host.getNumberOfPes();
		}
		// Divide the total resource usage by the number of hosts to get the average resource usage per host.
		double averageResourceUsage = totalResourceUsage / hostList.size();

		//System.out.println("Total usage: " + totalResourceUsage);

		// Calculate Fitness Function
		result = 1/makespan*averageResourceUsage;
		System.out.println("Fitness: " + result);
		return result;
	}
    
	private static Datacenter createDatacenter(String name){

		// Here are the steps needed to create a PowerDatacenter:
		// 1. We need to create a list to store
		//    our machine
		List<Host> hostList = new ArrayList<Host>();

		// 2. A Machine contains one or more PEs or CPUs/Cores.
		// In this example, it will have only one core.
		List<Pe> peList = new ArrayList<Pe>();

		int mips[] = {3500, 4000, 4500, 5000};

		// 3. Create PEs and add these into a list.
		peList.add(new Pe(0, new PeProvisionerSimple(mips[0]))); // need to store Pe id and MIPS Rating

		//4. Create Hosts with its id and list of PEs and add them to the list of machines
		int hostId=0;
		int ram = 50000; //host memory (MB)
		long storage = 1024000; //host storage
		int bw = 102400;

		hostList.add(
    			new Host(
    				hostId,
    				new RamProvisionerSimple(ram),
    				new BwProvisionerSimple(bw),
    				storage,
    				peList,
    				new VmSchedulerTimeShared(peList)
    			)
    	); // This is our first machine

		//create second machine in the Data center
		List<Pe> peList2 = new ArrayList<Pe>();

		peList2.add(new Pe(0, new PeProvisionerSimple(mips[1])));
		peList2.add(new Pe(1, new PeProvisionerSimple(mips[1])));

		hostId++;
		
		//specific
		ram = 100000; //host memory (MB)
		storage = 1024000; //host storage
		bw = 102400;

		hostList.add(
    			new Host(
    				hostId,
    				new RamProvisionerSimple(ram),
    				new BwProvisionerSimple(bw),
    				storage,
    				peList2,
    				new VmSchedulerTimeShared(peList2)
    			)
    	); // This is our second machine

		// Create third machine in the Data center
		List<Pe> peList3 = new ArrayList<Pe>();

		peList3.add(new Pe(0, new PeProvisionerSimple(mips[2])));
		peList3.add(new Pe(1, new PeProvisionerSimple(mips[2])));
		peList3.add(new Pe(2, new PeProvisionerSimple(mips[2])));

		hostId++;
		
		//specific
		ram = 150000; //host memory (MB)
		storage = 1024000; //host storage
		bw = 102400;

		hostList.add(
    			new Host(
    				hostId,
    				new RamProvisionerSimple(ram),
    				new BwProvisionerSimple(bw),
    				storage,
    				peList3,
    				new VmSchedulerTimeShared(peList3)
    			)
    	); // This is our third machine

		// Create our fourth machine in the Data center
		List<Pe> peList4 = new ArrayList<Pe>();

		peList4.add(new Pe(0, new PeProvisionerSimple(mips[3])));
		peList4.add(new Pe(1, new PeProvisionerSimple(mips[3])));
		peList4.add(new Pe(2, new PeProvisionerSimple(mips[3])));
		peList4.add(new Pe(3, new PeProvisionerSimple(mips[3])));

		hostId++;
		
		//specific
		ram = 200000; //host memory (MB)
		storage = 1024000; //host storage
		bw = 102400;

		hostList.add(
    			new Host(
    				hostId,
    				new RamProvisionerSimple(ram),
    				new BwProvisionerSimple(bw),
    				storage,
    				peList4,
    				new VmSchedulerTimeShared(peList4)
    			)
    	); // This is our fourth machine

		// 5. Create a DatacenterCharacteristics object that stores the
		//    properties of a data center: architecture, OS, list of
		//    Machines, allocation policy: time- or space-shared, time zone
		//    and its price (G$/Pe time unit).
		String arch = "x86";      // system architecture
		String os = "Linux";          // operating system
		String vmm = "Xen";
		double time_zone = 10.0;         // time zone this resource located
		double cost = 3.0;              // the cost of using processing in this resource
		double costPerMem = 0.05;		// the cost of using memory in this resource
		double costPerStorage = 0.001;	// the cost of using storage in this resource
		double costPerBw = 0.0;			// the cost of using bw in this resource
		LinkedList<Storage> storageList = new LinkedList<Storage>();	//we are not adding SAN devices by now

        DatacenterCharacteristics characteristics = new DatacenterCharacteristics(
                arch, os, vmm, hostList, time_zone, cost, costPerMem, costPerStorage, costPerBw);

		// 6. Finally, we need to create a PowerDatacenter object.
		Datacenter datacenter = null;
		try {
			datacenter = new Datacenter(name, characteristics, new VmAllocationPolicySimple(hostList), storageList, 0);
		} catch (Exception e) {
			e.printStackTrace();
		}

		return datacenter;
	}

	//We strongly encourage users to develop their own broker policies, to submit vms and cloudlets according
	//to the specific rules of the simulated scenario
	private static DatacenterBroker createBroker(){

		DatacenterBroker broker = null;
		try {
			broker = new DatacenterBroker("Broker");
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
		return broker;
	}
}
   
           
