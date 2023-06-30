#include <iostream>

// Function Prototypes

// Storyline functions. 
void beginning_choice() {
    std::cout << "As the Orks charge forward with reckless abandon, bellowing their war cries. The battalion holds fast, bolters primed. The air trembled with the roar of gunfire as the Orks closed in, you are faced with a choice.\n";
    std::cout << "What do you do?\n";
    std::cout << "1. Order your squad to leave the fortification and charge the Orks head on.\n";
    std::cout << "2. Order the squad to open fire.\n";
}

void leave_fortification() {
    std::cout << "You bellow the order of charge, your squad unquestioningly obeying your command bravely leap the machine gun emplacement and charge the Orks.\n"; 
    std::cout << "Captain Valerius, his storm bolter singing with righteous fury sees your squads courage and orders the rest of the fortification to charge. You have reached the first screaming Ork and faced with your second choice.\n";
    std::cout << "What do you do?\n";
    std::cout << "1. Unleash a flurry of blows with your power fist.\n";
    std::cout << "2. Release the fury of your storm bolter.\n";
}

void order_squad_open_fire() {
    std::cout << "The hail of bolter fire opens up on the advancing Orks. Yet, the Orks proved relentless, their brutish strength has allowed them to shrug off wounds that would cripple lesser beings.\n";
    std::cout << "The sheer numbers threaten to overwhelm your squad. Waves of green-skinned warriors surge forward, breaching your ranks. You are faced with a choice.\n";
    std::cout << "What do you do?\n";
    std::cout << "1. Order the squad to fall back.\n";
    std::cout << "2. Order the squad to hold the line and fight.\n";
}

void release_fury_of_storm_bolters() {
    std::cout << "You unleash a hail of bolter fire on the Ork. You spot a particularly large Ork, the Warboss, Grukka\n."; 
    std::cout << "You aim your storm bolter at the Warboss and unleash a flurry of bolts. The Warboss is hit, yet he continues to charge forward. You engage in hand to hand combat with the Warboss. You are faced with a choice.\n";
    std::cout << "What do you do?\n";
    std::cout << "1. Strike the Warboss with your power fist.\n";
    std::cout << "2. Order the squad to fall back.\n"; 
}

void power_fists() {
    std::cout << "The ground shakes beneath your feet as you engage in brutal hand to hand combat. The Ork swings his choppa with egregious force, aiming to cleave through your armor.\n";
    std::cout <<  "However your enhanced reflexes and training in close-quarters combat, deftly parry the Ork's attack. You strike his head with a thunderous blow, his skull shatters like glass. You are faced with a choice.\n";
    std::cout << "What do you do?\n";
    std::cout << "1. Continue engaging in hand to hand combat.\n";
    std::cout << "2. Order your flamer to unleash his mighty flamethrower.\n";
}

void fall_back() {
    std::cout << "You order the squad to fall back in an attempt to regroup. However, the Orks have breached your ranks and are in hot pursuit. Your team fails to regroup and are overwhelmed by the Orks.\n";
}

void hold_the_line() {
    std::cout << "Amidst the chaos, Librarian Tiberius, his mind aflame with psychic energy, channeled his powers to unleash a devastating psychic storm. Bolts of energy crackled through the air, incinerating Orks in their path.\n"; 
    std::cout << "The Warboss, Grukka, sensing the tide of battle turning against him, charged towards Tiberius with a thunderous roar. You are faced with a choice.\n";

    std::cout << "What do you do?\n";
    std::cout << "1. Alert Tiberius of the Warboss' charge.\n";
    std::cout << "2. Aim to strike the Warboss with your power fist.\n";
}

void losing_ending_1() {
    std::cout << "You continue to engage in hand to hand combat. The hoard relentless in their attack and you become overwhelmed. The Orks overcome your amour and you are torn apart by their choppas.\n";
}

void losing_ending_2() {
    std::cout << "Tiberius fails to hear your shout. The Warboss' charge eliminates the Librarian leaving your squad vulnerable to further attack, your squad becomes overwhelmed and perishes.\n";
}

void winning_ending() {
    std::cout << "Your order is the correct choice! The Orks are decimated, their corpses litter the ground, any foe remaining quickly loses moral and retreats. Your squad has held the flank. It is time to regroup with the rest of the battalion.\n";
}

//  Program Logic Functions. 
bool choice_checker(int choice){
    // Checks that the choice is valid. 
    if (choice <= 2 && choice > 0) {
        return true;
    }
    else {
        return false;
    }
}

bool end_checker(int *control_array_ptr, int choices) {
    // Checks to see if the player has made the maximum amount of choices possible.
    if (choices == 3){
        return true;
    }
    // Checks to see if the user has chosen the shorter path.
    else if (control_array_ptr[0] == 2 && control_array_ptr[1] == 1) {
        std::cout << "You have LOST.\n";
        return true;
    }
    // Else the user must not be finished yet.
    else {
        return false;
    }
}

int collect_user_input() {
    // Collects the input from the user. 

    int choice;
    std::cout << "Enter your choice: ";
    std::cin >> choice;
    // Checks that the user intput is valid.
    while (!choice_checker(choice)) {
        std::cout << "Invalid choice, please enter a valid choice: ";
        std::cin >> choice;
    }
    return choice;
}

void story_controller(int *control_array_ptr, int choices) {
    /* 
    Controls the flow of the story, takes an array of ints and the number of choices the user has made so far.
    Based on this information, the controller should call the correct function.
    */

    // Based on the number of choices, we dictate the output of the program using a switch statement.
    switch (choices) {
        // First choice.
        case 0:
            beginning_choice();
            break;

        // Second choice.
        case 1:
            if (control_array_ptr[0] == 1) {
                leave_fortification();
                break;
            }
            else if (control_array_ptr[0] == 2) {
                order_squad_open_fire();
                break;
            }
            else {
                std::cout << "Something went wrong, please restart the program.";
                break;
            }

        // Third choice.
        case 2: 
            if (control_array_ptr[0] == 1 && control_array_ptr[1] == 1) {
                power_fists();
                break;
            }
            else if (control_array_ptr[0] == 1 && control_array_ptr[1] == 2) {
                release_fury_of_storm_bolters();
                break;
            }
            else if (control_array_ptr[0] == 2 && control_array_ptr[1] == 1) {
                fall_back();
                break;
            }
            else if (control_array_ptr[0] == 2 && control_array_ptr[1] == 2) {
                hold_the_line();
                break;
            }
            else {
                std::cout << "Something went wrong, please restart the program.";
                break;
            }

        // Endings.
        case 3:
            // Most of the endings end positively, so check for bad endings else output a good ending.
            if (control_array_ptr[0] == 1 && control_array_ptr[1] == 1 && control_array_ptr[2] == 1) {
                losing_ending_1();
                break;
            }
            else if (control_array_ptr[0] == 1 && control_array_ptr[1] == 2 && control_array_ptr[2] == 2){
                fall_back();
                break;
            }
            else if (control_array_ptr[0] == 2 && control_array_ptr[1] == 2 && control_array_ptr[2] == 1){
                losing_ending_2();
                break;
            }
            else {
                winning_ending();
                break;
            }
        default:
            std::cout << "Something went wrong, please restart the program.";
            break;
    }
}

int main() {

    int user_decisions[3] = {0, 0, 0};
    int choices_made = 0;

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

    // First choice.
    story_controller(user_decisions, choices_made);
    user_decisions[choices_made] = collect_user_input();
    choices_made++;
    std::cout << "====================\n\n";

    // Second choice.
    story_controller(user_decisions, choices_made);
    user_decisions[choices_made] = collect_user_input();
    choices_made++;
    std::cout << "====================\n\n";

    // Third choice.
    story_controller(user_decisions, choices_made);
    if (end_checker(user_decisions, choices_made)) {
        return 0;
    }
    user_decisions[choices_made] = collect_user_input();
    choices_made++;
    std::cout << "====================\n\n";

    // Endings
    story_controller(user_decisions, choices_made);
    if (end_checker(user_decisions, choices_made)) {
        return 0;
    }
    else {
        return 1;
    }
}
