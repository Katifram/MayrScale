import torch
import torch.nn as nn

class NucleophilicityNet(nn.Module):
    def __init__(self):
        super(NucleophilicityNet, self).__init__()
        
        self.cnn_layers = nn.Sequential(
            nn.Conv1d(4, 64, 1),
            nn.ReLU(),
            nn.Conv1d(64, 128, 1),
            nn.ReLU(),
            nn.Conv1d(128, 256, 1),
            nn.ReLU(),
            nn.AdaptiveMaxPool1d(1)
        )

        self.fc_layers = nn.Sequential(
            nn.Linear(256 * 2 + 2, 128),
            nn.ReLU(),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 2)
        )

    def forward(self, schroedinger, energy_eigenvalues):

        batch_size = schroedinger.size(0)

        # Separate into homo and lumo
        schroedinger_lumo, schroedinger_homo = schroedinger[:,0,:,:], schroedinger[:,1,:,:]  

        schroedinger_lumo = schroedinger_lumo.transpose(1, 2) # (batch, 4, 1000)
        schroedinger_homo = schroedinger_homo.transpose(1, 2) # (batch, 4, 1000)

        schroedinger_lumo_cnn_output = self.cnn_layers(schroedinger_lumo).view(batch_size, -1)  # (batch, 256)
        schroedinger_homo_cnn_output = self.cnn_layers(schroedinger_homo).view(batch_size, -1)  # (batch, 256)

        combined_features = torch.cat([schroedinger_lumo_cnn_output, schroedinger_homo_cnn_output, energy_eigenvalues], dim=1)  # (batch, 514)

        output = self.fc_layers(combined_features) # (batch, 2)

        return output


# Example usage:
if __name__ == "__main__":

    model = NucleophilicityNet()

    schroedinger_solutions = torch.rand(1000, 2, 1000, 4)
    energy_eigenvalues = torch.rand(1000, 2)

    output = model(schroedinger_solutions, energy_eigenvalues)
    print(output.shape)
   
