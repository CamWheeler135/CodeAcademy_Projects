#include <iostream>

// Function Prototypes

bool choice_checker(int choice, int num_of_choices){
    if (choice <= num_of_choices && choice > 0) {
        return true;
    }
    else {
        return false;
    }
}

// Beginning Function
int beginning_choice() {
    // Text.
    std::cout << "As the Orks charge forward with reckless abandon, bellowing their war cries. The battalion holds fast, bolters primed. The air trembled with the roar of gunfire as the Orks closed in, you are faced with a choice.\n";
    // User Input
    int choice; 
    std::cout << "What do you choose?\n";
    std::cout << "1. Order your squad to leave the fortification and charge the Orks head on.\n" ;
    std::cout << "2. Order the squad to open fire.\n";
    std::cout << "Enter your choice: ";
    std::cin >> choice;
    // Check that input is acceptable.
    while (!choice_checker(choice, 2)) {
        std::cout << "Invalid choice, please enter a valid choice: ";
        std::cin >> choice;
    }
    return choice;
}

int first_layer_choice(int prev_layer_choice){
    if (prev_layer_choice == 1) {
        // Do this.
        return 0;
    }
    else if (prev_layer_choice == 2) {
        // Do that.
        return 0;
    }
    else {
        std::cout << "Error, choice passed into function does not have an action!\n";
        return 0;
    }
}

int main() {
    // Welcome message
    std::cout << "\nWelcome to my text based adventure game! My story is based on WarHammer 40K where the Space Marines face off in a battle against the Orks.\n";
    std::cout << "The game will offer you a series of choices, to select a choice simply enter the corresponding number when prompted and press enter!\n";
    std::cout << "I hope you enjoy!!!\n\n";

    // Story Begins.
    std::cout << "====================\n\n";
    std::cout << "The Story Begins!\n\n";
    std::cout << "In the grim darkness of the 41st millennium, the merciless forces of Chaos were not the only threat that plagued the Imperium of Man.\n"; 
    std::cout << "On a desolate world, a battalion of Space Marines, led by Captain Valerius of the Ultramarines, face an overwhelming horde of Orks led by the cunning Warboss Grukk.\n";
    std::cout << "You are a powerful Terminator in charge of the unit tasked with defending a flank of the frontline garrison against the impending ork attack.\n\n";
    std::cout << "====================\n\n";

    int first_choice = beginning_choice();

    return 0;
}